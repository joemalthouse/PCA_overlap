#!/usr/bin/env python3
"""
Pairwise Bhattacharyya coefficient heatmaps for comparing two
conformational ensembles in PCA space.

Discovers datasets from data/{label}_{state}_rep{N}.dat automatically.
Each cell (i, j) shows the 2D joint BC for (PC_i, PC_j); diagonal
cells show the 1D marginal BC.  Figures written to figures/.
"""

import glob
import re
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
FIGS = ROOT / "figures"

N_PCS     = 10
N_BINS_1D = 200
N_BINS_2D = 80

TITLES = {"cart": "Cartesian PCA", "dihed": "Dihedral PCA"}
STATE_LABELS = {"ox": "Oxidised", "red": "Reduced"}


def discover():
    datasets = {}
    for f in sorted(glob.glob(str(DATA / "*.dat"))):
        m = re.match(r"^(.+)_(.+)_rep(\d+)\.dat$", Path(f).name)
        if not m:
            continue
        label, state, _ = m.groups()
        datasets.setdefault(label, {}).setdefault(state, []).append(f)
    return datasets


def load(files):
    return np.vstack([np.loadtxt(f, comments="#")[:, 1:1+N_PCS] for f in files])


def bc_1d(x, y):
    edges = np.linspace(min(x.min(), y.min()),
                        max(x.max(), y.max()), N_BINS_1D + 1)
    px, _ = np.histogram(x, bins=edges, density=True)
    py, _ = np.histogram(y, bins=edges, density=True)
    return float(np.clip(np.sum(np.sqrt(px * py)) * (edges[1] - edges[0]), 0, 1))


def bc_2d(a, b, i, j):
    ei = np.linspace(min(a[:, i].min(), b[:, i].min()),
                     max(a[:, i].max(), b[:, i].max()), N_BINS_2D + 1)
    ej = np.linspace(min(a[:, j].min(), b[:, j].min()),
                     max(a[:, j].max(), b[:, j].max()), N_BINS_2D + 1)
    px, _, _ = np.histogram2d(a[:, i], a[:, j], bins=[ei, ej], density=True)
    py, _, _ = np.histogram2d(b[:, i], b[:, j], bins=[ei, ej], density=True)
    da = (ei[1] - ei[0]) * (ej[1] - ej[0])
    return float(np.clip(np.sum(np.sqrt(px * py)) * da, 0, 1))


def pairwise_matrix(a, b):
    n = min(N_PCS, a.shape[1], b.shape[1])
    mat = np.ones((n, n))
    for i in range(n):
        mat[i, i] = bc_1d(a[:, i], b[:, i])
        for j in range(i + 1, n):
            v = bc_2d(a, b, i, j)
            mat[i, j] = mat[j, i] = v
    return mat


def heatmap(mat, title, xlabel, ylabel, outfile):
    n = mat.shape[0]
    fig, ax = plt.subplots(figsize=(8, 7), constrained_layout=True)
    ax.spines["top"].set_visible(True)
    ax.spines["right"].set_visible(True)
    im = ax.imshow(mat, cmap="RdYlGn", vmin=0.5, vmax=1.0, origin="lower")
    for i in range(n):
        for j in range(n):
            v = mat[i, j]
            ax.text(j, i, f"{v:.2f}", ha="center", va="center",
                    fontsize=8, color="white" if v < 0.7 else "black")
    labels = [str(k + 1) for k in range(n)]
    ax.set_xticks(range(n))
    ax.set_xticklabels(labels)
    ax.set_yticks(range(n))
    ax.set_yticklabels(labels)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, shrink=0.8, label="Bhattacharyya coefficient")
    fig.savefig(outfile)
    plt.close(fig)


# ---- run ----

FIGS.mkdir(exist_ok=True)
datasets = discover()

for label, states in sorted(datasets.items()):
    names = sorted(states.keys())
    if len(names) != 2:
        continue
    state_a, state_b = names
    a, b = load(states[state_a]), load(states[state_b])
    mat = pairwise_matrix(a, b)

    pca_title = TITLES.get(label, f"{label} PCA")
    label_a = STATE_LABELS.get(state_a, state_a.capitalize())
    label_b = STATE_LABELS.get(state_b, state_b.capitalize())

    heatmap(mat,
            f"{pca_title} \u2014 {label_a} vs {label_b}",
            "Principal Component",
            "Principal Component",
            FIGS / f"bhattacharyya_{label}.pdf")

    print(f"{pca_title}: diagonal mean = {np.diag(mat).mean():.4f}, "
          f"off-diag min = {mat[np.triu_indices(mat.shape[0], 1)].min():.4f}")
