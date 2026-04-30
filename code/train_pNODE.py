import argparse
import os
import random
import time
from copy import deepcopy

import h5py
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from torchdiffeq import odeint
from torch.utils.data import Dataset, DataLoader


learning_rate = 0.001
num_epochs = 500
device = torch.device("cpu")
max_seed_attempts = 10

w_data = 0.5
w_ODE = 0.5

t3 = torch.tensor([0, 12, 24], dtype=torch.float32)
t5 = torch.tensor([0, 6, 12, 18, 24], dtype=torch.float32)
t10 = torch.tensor([0, 1, 2, 3, 7, 10, 12, 15, 18, 24], dtype=torch.float32)
t15 = torch.tensor([0, 1, 2, 3, 4, 6, 7, 9, 11, 13, 15, 17, 19, 21, 24], dtype=torch.float32)
t20 = torch.tensor([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24], dtype=torch.float32)
t25 = torch.tensor([0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24], dtype=torch.float32)

TIMEPOINTS_MAP = {
    "t3": t3,
    "t5": t5,
    "t10": t10,
    "t15": t15,
    "t20": t20,
    "t25": t25,
}


class CustomDataset(Dataset):
    def __init__(self, tensor):
        self.tensor = tensor

    def __len__(self):
        return len(self.tensor)

    def __getitem__(self, idx):
        input_data = self.tensor[idx, 0, :]
        output_data = self.tensor[idx, :, :]
        return input_data, output_data


def loss_fn(pred_y, y, pred_y_ode, dydt_pred, model):
    data_loss = torch.mean((y - pred_y) ** 2)
    glv_rhs = pred_y_ode * (model.alpha + torch.matmul(pred_y_ode, model.beta.T))
    dydt_residual = dydt_pred - glv_rhs
    ode_loss = torch.mean(dydt_residual ** 2)
    total_loss = w_data * data_loss + w_ODE * ode_loss
    return total_loss


class ODEFunc(nn.Module):
    def __init__(self, p):
        super().__init__()
        self.P = p

        self.alpha = nn.Parameter(torch.rand(p, device=device))
        self.beta = nn.Parameter(self._generate_A(p))
        self.net = nn.Sequential(
            nn.Linear(p, 2 * p),
            nn.Tanh(),
            nn.Linear(2 * p, p),
        )

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, mean=0.0, std=0.1)
                nn.init.constant_(m.bias, val=0.0)

    def forward(self, t, y):
        I = (y > 0).int()
        return torch.mul(self.net(y), I)

    def _generate_A(self, p):
        A_data = torch.normal(mean=0.0, std=0.1, size=(p, p), device=device)
        for i in range(p):
            A_data[i, i] = -1.0
        return A_data

    def constrained(self):
        with torch.no_grad():
            self.alpha.data = torch.relu(self.alpha.data)


def set_seed(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def parse_args():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=10)
    parser.add_argument("--s_label", type=str, required=True, choices=["s10", "s20"])
    parser.add_argument("--t_label", type=str, required=True, choices=["t5", "t10", "t15", "t20", "t25"])
    parser.add_argument("--input_root", type=str, required=True)
    parser.add_argument("--output_root", type=str, required=True)
    parser.add_argument("--start_iter", type=int, default=1)
    parser.add_argument("--end_iter", type=int, default=10)
    return parser.parse_args()


def main():
    args = parse_args()

    P = args.p
    s = args.s_label
    t = args.t_label

    timepoints = TIMEPOINTS_MAP[t].to(device)
    t_ode = torch.arange(0, 25, dtype=torch.float32, device=device)

    for iter_num in range(args.start_iter, args.end_iter + 1):
        current_seed = 1
        successful_train = False

        while (not successful_train) and (current_seed <= max_seed_attempts):
            try:
                set_seed(current_seed)

                input_dir = os.path.join(
                    args.input_root,
                    f"p{P}",
                    s,
                    t,
                    f"iter{iter_num}",
                )
                output_dir = os.path.join(
                    args.output_root,
                    f"p{P}",
                    "training",
                    s,
                    t,
                    f"iter{iter_num}",
                )

                os.makedirs(output_dir, exist_ok=True)

                data_file = os.path.join(input_dir, "true_initial_state.h5")
                with h5py.File(data_file, "r") as hdf:
                    data = hdf["true_initial_state"][:]
                    tensor_data = torch.tensor(data, dtype=torch.float32).permute(0, 2, 1)

                dataset = CustomDataset(tensor_data)
                batch_s = int(int(s[1:]) * 0.2)
                dataloader = DataLoader(dataset, batch_size=batch_s, shuffle=False)

                model = ODEFunc(P).to(device)
                optimizer = optim.Adam(model.parameters(), lr=learning_rate)

                best_loss = float("inf")
                best_model = None
                loss_history = []

                avg_loss_file_path = os.path.join(output_dir, "average_loss.txt")
                with open(avg_loss_file_path, "w") as f:
                    f.write("")

                start_time = time.time()

                for epoch in range(num_epochs):
                    for input_data, output_data in dataloader:
                        input_data = input_data.to(device).float()
                        output_data = output_data.to(device)

                        optimizer.zero_grad()

                        pred_y = odeint(model, input_data, timepoints).transpose(1, 0)
                        pred_y_ode = odeint(model, input_data, t_ode).transpose(1, 0)
                        dydt_pred = model(0, pred_y_ode)

                        loss = loss_fn(pred_y, output_data, pred_y_ode, dydt_pred, model)
                        loss.backward()
                        optimizer.step()
                        model.constrained()

                        if loss.item() < best_loss:
                            best_loss = deepcopy(loss.item())
                            best_model = deepcopy(model.state_dict())

                        loss_history.append(loss.item())

                    avg_loss = np.mean(loss_history[-len(dataloader):])
                    print(f"[{s}, {t}, iter{iter_num}] Epoch [{epoch + 1}/{num_epochs}], Average Loss: {avg_loss:.4f}")
                    with open(avg_loss_file_path, "a") as f:
                        f.write(f"Epoch [{epoch + 1}/{num_epochs}], Average Loss: {avg_loss:.4f}\n")

                total_time = time.time() - start_time

                torch.save(best_model, os.path.join(output_dir, "best_model.pth"))
                np.savetxt(os.path.join(output_dir, "loss_history.txt"), loss_history)

                runtime_file_path = os.path.join(output_dir, "runtime.txt")
                with open(runtime_file_path, "w") as f:
                    f.write(f"Model: {s}/{t}/iter{iter_num}. Total runtime: {total_time:.2f} seconds\n")

                print(f"Model trained for {s}/{t}/iter{iter_num}. Total runtime: {total_time:.2f} seconds")
                successful_train = True
                print(f"Iter {iter_num} trained successfully with seed {current_seed}.")

            except AssertionError as e:
                if "underflow in dt" in str(e):
                    print(f"Underflow error on Iter {iter_num} with seed {current_seed}. Trying a different seed...")
                    current_seed += 1
                else:
                    raise

        if not successful_train:
            print(f"Failed to train Iter {iter_num} after trying {max_seed_attempts} seeds.")
            break


if __name__ == "__main__":
    main()