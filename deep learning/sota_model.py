import torch
import torch.nn as nn

#  Define model with 7 hidden layers
class DeepFFNN(nn.Module):
    def __init__(self):
        super(DeepFFNN, self).__init__()
        self.sigmoid = nn.Sigmoid()
        self.layers = nn.Sequential(
            nn.Linear(29,29),
            nn.Sigmoid(),
            nn.Linear(29,29),
            nn.Sigmoid(),
            nn.Linear(29,29),
            nn.Sigmoid(),
            nn.Linear(29,29),
            nn.Sigmoid(),
            nn.Linear(29,29),
            nn.Sigmoid(),
            nn.Linear(29,29),
            nn.Sigmoid(),
            nn.Linear(29,29),
            nn.Sigmoid(),
            nn.Linear(29, 1)
        )

    def forward(self, x):
        return self.layers(x)
