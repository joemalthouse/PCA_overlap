#!/usr/bin/env python3
"""
KDE contour plots of PC1-3 pairwise projections for each state and
PCA type.

Discovers datasets from data/{label}_{state}_rep{N}.dat automatically.
Variance fractions are read from data/{label}_eigenfrac.dat if present.

Figures are written to figures/.
"""

import glob
import re
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from scipy.stats import gaussian_kde
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]
DATA = ROOT / "data"
FIGS = ROOT / "figures"

SUBSAMPLE = 10
GRID_PTS  = 200
LEVELS    = np.arange(0.95, 0.04, -0.05)

STATE_STYLE = {
    "ox":  {"cmap": plt.cm.Blues, "label": "Oxidised"},
    "red": {"cmap": plt.cm.Reds,  "label": "Reduced"},
}
DEFAULT_CMAPS = [plt.cm.Blues, plt.cm.Reds, plt.cm.Greens, plt.cm.Purples]

LABEL_TITLES = {"cart": "Cartesian PCA", "dihed": "Dihedral PCA"}


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
    return np.vstack([np.loadtxt(f, comments="#")[:, 1:4] for f in files])


def load_eigenfrac(label):
    path = DATA / f"{label}_eigenfrac.dat"
    if not path.exists():
        return None
    return np.loadtxt(path, comments="#")[:, 1]


def density_thresholds(Z):
    flat = Z.ravel()
    idx = np.argsort(flat)[::-1]
    cum = np.cumsum(flat[idx])
    cum /= cum[-1]
    return np.array([flat[idx[min(np.searchsorted(cum, lev), len(flat) - 1)]]
                     for lev in LEVELS])


def make_overlay(all_sub, grid_lo, grid_hi, varfrac, state_info, title):
    """Overlay KDE contours for all states on the same axes."""
    pairs = [(0, 1), (0, 2), (1, 2)]
    if varfrac is not None:
        pc_labels = [f"PC{i+1} ({varfrac[i]*100:.1f}%)" for i in range(3)]
    else:
        pc_labels = [f"PC{i+1}" for i in range(3)]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

    for ax, (i, j) in zip(axes, pairs):
        xi = np.linspace(grid_lo, grid_hi, GRID_PTS)
        xj = np.linspace(grid_lo, grid_hi, GRID_PTS)
        XI, XJ = np.meshgrid(xi, xj)
        grid_coords = np.vstack([XI.ravel(), XJ.ravel()])

        for data_sub, cmap, slabel in state_info:
            kde = gaussian_kde(data_sub[:, [i, j]].T)
            Z = kde(grid_coords).reshape(XI.shape)

            thresh = np.sort(density_thresholds(Z))
            n = len(thresh)
            colors = [cmap(0.3 + 0.6 * k / (n - 1)) for k in range(n)]

            ax.contour(XI, XJ, Z, levels=thresh,
                       colors=colors, linewidths=1.2, alpha=0.75)
            ax.contourf(XI, XJ, Z, levels=[thresh[0], Z.max()],
                        colors=[cmap(0.15)], alpha=0.10)

        ax.set_xlabel(pc_labels[i], fontsize=24)
        ax.set_ylabel(pc_labels[j], fontsize=24)
        ax.tick_params(labelsize=18)
        ax.set_facecolor("#f7f7f7")

    from matplotlib.lines import Line2D
    handles = [Line2D([0], [0], color=cmap(0.7), lw=2, label=slabel)
               for _, cmap, slabel in state_info]
    axes[2].legend(handles=handles, loc="upper right", fontsize=18,
                   framealpha=0.9)

    fig.suptitle(title, fontsize=28)
    fig.subplots_adjust(left=0.06, right=0.97, bottom=0.12, top=0.90,
                        wspace=0.30)
    return fig


def make_figure(data_sub, grid_lo, grid_hi, varfrac, cmap, title):
    pairs = [(0, 1), (0, 2), (1, 2)]
    if varfrac is not None:
        pc_labels = [f"PC{i+1} ({varfrac[i]*100:.1f}%)" for i in range(3)]
    else:
        pc_labels = [f"PC{i+1}" for i in range(3)]

    fig, axes = plt.subplots(1, 3, figsize=(17, 5.2))

    for ax, (i, j) in zip(axes, pairs):
        xi = np.linspace(grid_lo, grid_hi, GRID_PTS)
        xj = np.linspace(grid_lo, grid_hi, GRID_PTS)
        XI, XJ = np.meshgrid(xi, xj)

        kde = gaussian_kde(data_sub[:, [i, j]].T)
        Z = kde(np.vstack([XI.ravel(), XJ.ravel()])).reshape(XI.shape)

        thresh = np.sort(density_thresholds(Z))
        n = len(thresh)
        colors = [cmap(0.3 + 0.6 * k / (n - 1)) for k in range(n)]

        ax.contour(XI, XJ, Z, levels=thresh,
                   colors=colors, linewidths=1.2, alpha=0.85)
        ax.contourf(XI, XJ, Z, levels=[thresh[0], Z.max()],
                    colors=[cmap(0.15)], alpha=0.15)
        ax.set_xlabel(pc_labels[i], fontsize=24)
        ax.set_ylabel(pc_labels[j], fontsize=24)
        ax.tick_params(labelsize=18)
        ax.set_facecolor("#f7f7f7")

    fig.suptitle(title, fontsize=28)
    fig.subplots_adjust(left=0.06, right=0.97, bottom=0.12, top=0.90,
                        wspace=0.30)
    return fig


# ---- run ----

FIGS.mkdir(exist_ok=True)
datasets = discover()
if not datasets:
    print(f"No data found in {DATA}/")
    raise SystemExit(1)

for label, states in sorted(datasets.items()):
    varfrac = load_eigenfrac(label)
    all_data = {s: load(files) for s, files in states.items()}

    pooled = np.vstack(list(all_data.values()))[::SUBSAMPLE]
    lo, hi = pooled.min(), pooled.max()
    pad = 0.1 * (hi - lo)
    grid_lo, grid_hi = lo - pad, hi + pad

    pca_title = LABEL_TITLES.get(label, f"{label} PCA")

    for idx, (state, raw) in enumerate(sorted(all_data.items())):
        if state in STATE_STYLE:
            cmap = STATE_STYLE[state]["cmap"]
            state_label = STATE_STYLE[state]["label"]
        else:
            cmap = DEFAULT_CMAPS[idx % len(DEFAULT_CMAPS)]
            state_label = state.capitalize()

        title = f"{pca_title} \u2014 {state_label}"
        fig = make_figure(raw[::SUBSAMPLE], grid_lo, grid_hi, varfrac,
                          cmap, title)
        out = FIGS / f"kde_{label}_{state}.pdf"
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        print(out)

    # Overlay figure with both states
    if len(all_data) == 2:
        state_info = []
        for sidx, (state, raw) in enumerate(sorted(all_data.items())):
            if state in STATE_STYLE:
                cmap = STATE_STYLE[state]["cmap"]
                slabel = STATE_STYLE[state]["label"]
            else:
                cmap = DEFAULT_CMAPS[sidx % len(DEFAULT_CMAPS)]
                slabel = state.capitalize()
            state_info.append((raw[::SUBSAMPLE], cmap, slabel))

        fig = make_overlay(None, grid_lo, grid_hi, varfrac, state_info,
                           f"{pca_title} \u2014 Overlay")
        out = FIGS / f"kde_{label}_overlay.pdf"
        fig.savefig(out, bbox_inches="tight", pad_inches=0.05)
        plt.close(fig)
        print(out)
