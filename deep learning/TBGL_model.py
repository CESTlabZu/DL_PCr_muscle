
"""
Two-branch Global-Local Model
"""

from __future__ import annotations

import math
import torch
import torch.nn as nn
import torch.nn.functional as F

class ConvDenoiser(nn.Module):
    """Simple 1D Conv denoiser.

    Parameters
    ----------
    in_channels : int, default=2
        Number of input channels. Output has the same channel count.
    """

    def __init__(self, in_channels: int = 2) -> None:
        super().__init__()
        self.net = nn.Sequential(
            nn.Conv1d(in_channels, 128, kernel_size=5, padding=2),
            nn.ELU(),
            nn.Conv1d(128, 128, kernel_size=3, padding=1),
            nn.ELU(),
            nn.Conv1d(128, 128, kernel_size=3, padding=1),
            nn.ELU(),
            # Output remains `in_channels` (denoised signal)
            nn.Conv1d(128, in_channels, kernel_size=3, padding=1),
        )

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Apply the convolutional denoiser.

        Parameters
        ----------
        x : torch.Tensor
            Shape: (B, C=in_channels, L)
        """
        return self.net(x)


class SecondDerivativeFiniteDiff(nn.Module):
    """Second derivative (finite-difference) layer.

    Computes d^2 x / d n^2 along the last dimension via a fixed kernel.
    """

    def __init__(self) -> None:
        super().__init__()
        kernel = torch.tensor([[-1.0, 2.0, -1.0]], dtype=torch.float32)  # (1, 3)
        kernel = kernel.unsqueeze(0)  # shape: (1, 1, 3)
        self.register_buffer("kernel", kernel)

    def forward(self, x: torch.Tensor) -> torch.Tensor:
        """Apply second-derivative filter channel-wise.

        Parameters
        ----------
        x : torch.Tensor
            Shape: (B, C, L)
        """
        b, c, l = x.shape
        x = x.view(b * c, 1, l)
        d2x = F.conv1d(x, self.kernel, padding=1)
        return d2x.view(b, c, l)


class TwoBranchRegressor(nn.Module):
    """Two-branch model: local + global streams with denoising front-ends.

    Outputs two scalars (fs, ksw) and returns denoised inputs

    Parameters
    ----------
    local_feat_dim : int, default=128
    global_feat_dim : int, default=128
    fusion_dim : int, default=128
    """

    def __init__(
        self,
        local_feat_dim: int = 128,
        global_feat_dim: int = 128,
        fusion_dim: int = 128,
    ) -> None:
        super().__init__()

        # Denoisers (2-channel input)
        self.global_denoiser = ConvDenoiser(in_channels=2)
        self.local_denoiser = ConvDenoiser(in_channels=2)

        # Finite-difference feature (applied to local branch)
        self.diff_layer = SecondDerivativeFiniteDiff()

        # Local feature extractor
        self.local_encoder = nn.Sequential(
            nn.Conv1d(2, 32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(64, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
        )
        self.local_fc = nn.Linear(64, local_feat_dim)

        # Global feature extractor
        self.global_encoder = nn.Sequential(
            nn.Conv1d(2, 32, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(32, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.Conv1d(64, 64, kernel_size=3, padding=1),
            nn.ReLU(),
            nn.AdaptiveAvgPool1d(1),
        )
        self.global_fc = nn.Linear(64, global_feat_dim)

        # Fusion and regression heads
        self.fusion_fc = nn.Sequential(
            nn.Linear(local_feat_dim + global_feat_dim, fusion_dim),
            nn.ReLU(),
            nn.Linear(fusion_dim, 64),
            nn.ReLU(),
        )
        self.fs_head = nn.Linear(64, 1)
        self.ksw_head = nn.Linear(64, 1)

    def forward(
        self, global_input: torch.Tensor, local_input: torch.Tensor
    ) -> tuple[torch.Tensor, torch.Tensor, torch.Tensor, torch.Tensor]:
        """Forward pass.

        Parameters
        ----------
        global_input : torch.Tensor
            Shape: (B, 2, Lg)
        local_input : torch.Tensor
            Shape: (B, 2, Ll)

        Returns
        -------
        fs_reg : torch.Tensor
            Shape: (B, 1)
        ksw_reg : torch.Tensor
            Shape: (B, 1)
        denoised_local : torch.Tensor
            Denoised local input, shape: (B, 2, Ll)
        denoised_global : torch.Tensor
            Denoised global input, shape: (B, 2, Lg)
        """
        # Denoise inputs
        denoised_global = self.global_denoiser(global_input)
        denoised_local = self.local_denoiser(local_input)

        # Second-derivative feature for local branch
        derived_local = self.diff_layer(denoised_local)

        # Local pathway
        local_feat = self.local_encoder(derived_local).squeeze(-1)
        local_feat = self.local_fc(local_feat)

        # Global pathway
        global_feat = self.global_encoder(denoised_global).squeeze(-1)
        global_feat = self.global_fc(global_feat)

        # Fuse and regress
        fused_feat = torch.cat([local_feat, global_feat], dim=1)
        hidden = self.fusion_fc(fused_feat)
        fs_reg = self.fs_head(hidden)
        ksw_reg = self.ksw_head(hidden)

        return fs_reg, ksw_reg, denoised_local, denoised_global
