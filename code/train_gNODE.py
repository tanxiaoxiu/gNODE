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


def loss_fn(pred_y, y):
    return torch.mean((y - pred_y) ** 2)


class ODEFunc(nn.Module):
    def __init__(self, p):
        super().__init__()
        self.P = p
        self.alpha = nn.Parameter(torch.rand(p, device=device))
        self.beta = nn.Parameter(self._generate_A(p))

    def forward(self, t, ys):
        dydts_list = []
        for y in ys:
            dydt = y * (self.alpha + torch.matmul(self.beta, y))
            dydt = dydt * (y > 0).float()
            dydts_list.append(dydt)
        dydts = torch.stack(dydts_list)
        return dydts

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

    for iter_num in range(args.start_iter, args.end_iter + 1):
        current_seed = 1
        successful_train = False

        while not successful_train and current_seed < max_seed_attempts:
            try:
                set_seed(current_seed)

                input_file_path = os.path.join(
                    args.input_root,
                    f"p{P}",
                    s,
                    t,
                    f"iter{iter_num}",
                    "true_initial_state.h5",
                )

                output_dir_path = os.path.join(
                    args.output_root,
                    f"p{P}",
                    "training",
                    s,
                    t,
                    f"iter{iter_num}",
                )

                os.makedirs(output_dir_path, exist_ok=True)

                with h5py.File(input_file_path, "r") as hdf:
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

                avg_loss_file_path = os.path.join(output_dir_path, "average_loss.txt")
                with open(avg_loss_file_path, "w") as file:
                    file.write("")

                start_time = time.time()

                for epoch in range(num_epochs):
                    for batch_idx, (input_data, output_data) in enumerate(dataloader):
                        input_data, output_data = input_data.to(device), output_data.to(device)
                        optimizer.zero_grad()
                        input_data = input_data.float()

                        pred_y = odeint(model, input_data, timepoints)
                        pred_y = pred_y.transpose(1, 0)
                        loss = loss_fn(pred_y, output_data)

                        loss.backward()
                        optimizer.step()
                        model.constrained()

                        if loss.item() < best_loss:
                            best_loss = deepcopy(loss.item())
                            best_model = deepcopy(model.state_dict())

                        loss_history.append(loss.item())

                    avg_loss = np.mean(loss_history[-len(dataloader):])
                    print(f"[{s}, {t}, iter{iter_num}] Epoch [{epoch + 1}/{num_epochs}], Average Loss: {avg_loss:.4f}")
                    with open(avg_loss_file_path, "a") as file:
                        file.write(f"Epoch [{epoch + 1}/{num_epochs}], Average Loss: {avg_loss:.4f}\n")

                end_time = time.time()
                total_time = end_time - start_time

                torch.save(best_model, os.path.join(output_dir_path, "best_model.pth"))
                np.savetxt(os.path.join(output_dir_path, "loss_history.txt"), loss_history)

                runtime_file_path = os.path.join(output_dir_path, "runtime.txt")
                with open(runtime_file_path, "w") as file:
                    file.write(f"Model: {s}/{t}/iter{iter_num}. Total runtime: {total_time:.2f} seconds\n")

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