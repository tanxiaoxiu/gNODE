import argparse
import os
import random
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy

import h5py
import numpy as np
import torch
import torch.nn as nn
import torch.optim as optim
from sklearn.model_selection import KFold
from torch.utils.data import DataLoader, Dataset
from torchdiffeq import odeint

learning_rate = 0.001
num_epochs = 500
device = torch.device("cpu")

t3 = torch.tensor([0, 12, 24], dtype=torch.float32)
t5 = torch.tensor([0, 6, 12, 18, 24], dtype=torch.float32)
t10 = torch.tensor([0, 1, 2, 3, 7, 10, 12, 15, 18, 24], dtype=torch.float32)
t15 = torch.tensor([0, 1, 2, 3, 4, 6, 7, 9, 11, 13, 15, 17, 19, 21, 24], dtype=torch.float32)
t20 = torch.tensor(
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 16, 18, 20, 22, 24],
    dtype=torch.float32,
)
t25 = torch.tensor(
    [0, 1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24],
    dtype=torch.float32,
)

timepoint_map = {
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
    def __init__(self, p_dim, run_device):
        super().__init__()
        self.p_dim = p_dim
        self.net = nn.Sequential(
            nn.Linear(p_dim, 2 * p_dim),
            nn.Tanh(),
            nn.Linear(2 * p_dim, p_dim),
        )

        for m in self.net.modules():
            if isinstance(m, nn.Linear):
                nn.init.normal_(m.weight, mean=0.0, std=0.1)
                nn.init.constant_(m.bias, val=0.0)

    def forward(self, t, y):
        indicator = (y > 0).int()
        return torch.mul(self.net(y), indicator)


def set_seed(seed):
    torch.manual_seed(seed)
    if torch.cuda.is_available():
        torch.cuda.manual_seed_all(seed)
    np.random.seed(seed)
    random.seed(seed)
    torch.backends.cudnn.deterministic = True
    torch.backends.cudnn.benchmark = False


def evaluate_full_validation_set(model, val_dataset, timepoints, run_device, output_dir, fold_idx):
    model.eval()
    predictions = []
    val_targets = []

    with torch.no_grad():
        for i in range(len(val_dataset)):
            val_input, val_output = val_dataset[i]
            val_input = val_input.to(run_device).float()
            val_output = val_output.to(run_device)

            pred_val_y = odeint(model, val_input.unsqueeze(0), timepoints)[:, 0, :]  # [T, p]

            predictions.append(pred_val_y.cpu().numpy())
            val_targets.append(val_output.cpu().numpy())

            pred_val_y_file_path = os.path.join(
                output_dir, f"pred_val_y_fold_{fold_idx}_individual_{i}.txt"
            )
            np.savetxt(pred_val_y_file_path, pred_val_y.cpu().numpy(), fmt="%f")

    predictions = np.stack(predictions, axis=0)
    val_targets = np.stack(val_targets, axis=0)

    if predictions.shape != val_targets.shape:
        raise ValueError(
            f"Shape mismatch in validation: predictions {predictions.shape} vs targets {val_targets.shape}"
        )

    rmse = np.sqrt(np.mean((val_targets - predictions) ** 2))
    relative_rmse = rmse / np.sqrt(np.mean(val_targets ** 2))
    return relative_rmse


def train_and_evaluate(train_data, val_data, run_device, batch_size, output_dir, fold_idx, timepoints, p_dim):
    set_seed(1)

    if not isinstance(train_data, torch.Tensor):
        train_data = torch.tensor(train_data, dtype=torch.float32)
    if not isinstance(val_data, torch.Tensor):
        val_data = torch.tensor(val_data, dtype=torch.float32)

    train_dataset = CustomDataset(train_data)
    val_dataset = CustomDataset(val_data)
    train_loader = DataLoader(train_dataset, batch_size=batch_size, shuffle=True)

    model = ODEFunc(p_dim, run_device).to(run_device)
    optimizer = optim.Adam(model.parameters(), lr=learning_rate)

    best_model = None
    best_loss = float("inf")
    loss_history = []

    avg_loss_file = os.path.join(output_dir, f"average_loss_fold_{fold_idx}.txt")
    with open(avg_loss_file, "w", encoding="utf-8") as file:
        file.write("")

    for epoch in range(num_epochs):
        epoch_losses = []

        for input_data, output_data in train_loader:
            input_data = input_data.to(run_device).float()
            output_data = output_data.to(run_device)

            optimizer.zero_grad()
            pred_y = odeint(model, input_data, timepoints)
            pred_y = pred_y.transpose(1, 0)
            loss = loss_fn(pred_y, output_data)

            loss.backward()
            optimizer.step()

            epoch_losses.append(loss.item())

            if loss.item() < best_loss:
                best_loss = deepcopy(loss.item())
                best_model = deepcopy(model.state_dict())

        avg_loss = np.mean(epoch_losses)
        with open(avg_loss_file, "a", encoding="utf-8") as file:
            file.write(f"Epoch [{epoch + 1}/{num_epochs}], Average Loss: {avg_loss:.4f}\n")

        loss_history.extend(epoch_losses)

    if best_model is None:
        raise RuntimeError(f"Fold {fold_idx}: best_model was never set.")

    model_path = os.path.join(output_dir, f"best_model_fold_{fold_idx}.pth")
    torch.save(best_model, model_path)

    loss_file_path = os.path.join(output_dir, f"loss_history_fold_{fold_idx}.txt")
    np.savetxt(loss_file_path, loss_history, fmt="%f")

    model.load_state_dict(best_model)
    relative_rmse = evaluate_full_validation_set(
        model, val_dataset, timepoints, run_device, output_dir, fold_idx
    )

    print(f"Fold {fold_idx} trained successfully.")
    return fold_idx, relative_rmse


def k_fold_cross_validation(data, num_folds, run_device, batch_size, output_dir, timepoints, p_dim):
    kf = KFold(n_splits=num_folds, shuffle=True, random_state=1)
    rmses = [None] * num_folds

    with ProcessPoolExecutor(max_workers=5) as executor:
        futures = []
        for fold_idx, (train_index, val_index) in enumerate(kf.split(data)):
            train_data = data[train_index]
            val_data = data[val_index]
            futures.append(
                executor.submit(
                    train_and_evaluate,
                    train_data,
                    val_data,
                    run_device,
                    batch_size,
                    output_dir,
                    fold_idx,
                    timepoints,
                    p_dim,
                )
            )

        for future in futures:
            fold_idx, relative_rmse = future.result()
            rmses[fold_idx] = relative_rmse

    rmses_file_path = os.path.join(output_dir, "k_fold_rrmses.txt")
    np.savetxt(rmses_file_path, rmses, fmt="%f")
    avg_rmse = np.mean(rmses)
    return avg_rmse


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--p", type=int, default=10)
    parser.add_argument("--sparsity_label", type=str, required=True)
    parser.add_argument("--s_label", type=str, required=True)
    parser.add_argument("--t_label", type=str, required=True)
    parser.add_argument("--input_root", type=str, required=True)
    parser.add_argument("--output_root", type=str, required=True)
    parser.add_argument("--iter_start", type=int, default=1)
    parser.add_argument("--iter_end", type=int, default=10)
    args = parser.parse_args()

    p_dim = args.p
    if args.t_label not in timepoint_map:
        raise ValueError(f"Unsupported t_label: {args.t_label}")

    timepoints = timepoint_map[args.t_label].to(device)

    for iter_num in range(args.iter_start, args.iter_end + 1):
        input_file_path = os.path.join(
            args.input_root, f"p{p_dim}", args.s_label, args.t_label, f"iter{iter_num}"
        )

        output_dir_path = os.path.join(
            args.output_root,
            f"p{p_dim}",
            "training",
            args.s_label,
            args.t_label,
            f"iter{iter_num}",
        )
        os.makedirs(output_dir_path, exist_ok=True)

        data_file = os.path.join(input_file_path, "true_initial_state.h5")
        if not os.path.exists(data_file):
            raise FileNotFoundError(f"Missing input file: {data_file}")

        with h5py.File(data_file, "r") as hdf:
            data = hdf["true_initial_state"][:]
            tensor_data = torch.tensor(data, dtype=torch.float32).permute(0, 2, 1)

        num_folds = 5
        batch_size = max(1, int(tensor_data.shape[0] * 0.2))

        avg_rmse = k_fold_cross_validation(
            tensor_data,
            num_folds,
            device,
            batch_size,
            output_dir_path,
            timepoints,
            p_dim,
        )

        print(
            f"sparsity={args.sparsity_label}, subject={args.s_label}, "
            f"time={args.t_label}, iter={iter_num}, "
            f"Average Relative RMSE: {avg_rmse:.4f}"
        )

        avg_rmse_file_path = os.path.join(output_dir_path, "average_relative_rmse.txt")
        with open(avg_rmse_file_path, "w", encoding="utf-8") as file:
            file.write(f"Average Relative RMSE: {avg_rmse:.4f}\n")


if __name__ == "__main__":
    main()