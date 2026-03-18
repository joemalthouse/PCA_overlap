# PCA Ensemble Comparison

Compares oxidised and reduced m4D2 ensembles in PCA space using two approaches:

1. **Bhattacharyya coefficient heatmaps** — pairwise overlap between principal
   components, with 2D joint BC off-diagonal and 1D marginal on the diagonal.
2. **KDE contour plots** — probability density in the PC1–3 plane for each
   state individually and as an overlay.

Both cartesian (Cα) and dihedral (φ/ψ) PCA are included.

## Folder structure

```
PCA/
├── README.md
├── run.sh                            # Sequential fallback (no SLURM)
├── scripts/
│   ├── cpptraj/
│   │   ├── 02_cartesian_pca.cpptraj  # Cα covariance → eigendecomp → projection
│   │   └── 03_dihedral_pca.cpptraj   # φ/ψ covariance → eigendecomp → projection
│   └── python/
│       ├── bhattacharyya.py          # BC heatmaps
│       └── kde_contours.py           # KDE contour plots
├── slurm/
│   ├── 02_cart_pca.slurm
│   ├── 03_dihed_pca.slurm
│   ├── 04_figures.slurm
│   └── submit_all.sh                # Dependency-chained submission
├── data/                             # Projections and eigenfractions (generated)
├── figures/                          # PDFs (generated)
└── logs/                             # SLURM logs
```

## Prerequisites

The cpptraj scripts expect:

- **Amber topologies** (`.prmtop`) for each state
- **Cα-aligned trajectory streams** (NetCDF) for each replica of each state
- **Joint Cα average structure** (PDB) — used as the RMS reference for
  cartesian PCA (not needed for dihedral PCA)

Paths to these files are set at the top of each `.cpptraj` script.

## Running

**With SLURM** (cartesian and dihedral PCA run in parallel):

```bash
bash PCA/slurm/submit_all.sh
```

**Without SLURM** (sequential):

```bash
bash PCA/run.sh
```

Both clean previous outputs before starting.

## Pipeline

1. **Cartesian PCA** — loads all 20 aligned streams (10 ox + 10 red), fits to
   the joint Cα average, builds the Cα covariance matrix, diagonalises for 20
   eigenvectors, then projects each replica individually.
2. **Dihedral PCA** — same streams, computes φ/ψ via `multidihedral`, builds
   the dihedral covariance matrix (`matrix dihcovar`), diagonalises, then
   projects each replica's dihedrals individually.
3. **Figures** — Python scripts auto-discover whatever `.dat` files are in
   `data/` and produce the heatmaps and contour plots.

## Data format

Projection files follow the naming convention `{label}_{state}_rep{N}.dat`:

```
#Frame   Mode1    Mode2    Mode3   ...
     1   1.424   -0.121   -5.235  ...
     2   0.397   -1.407   -5.780  ...
```

Column 0 is the frame index; columns 1+ are PC projections. Eigenfraction
files (`{label}_eigenfrac.dat`) contain variance explained per mode and are
needed for axis labels on the KDE plots.

## Configuration

| Variable | Default | Script | Purpose |
|----------|---------|--------|---------|
| `N_PCS` | 10 | bhattacharyya.py | PCs included in BC analysis |
| `N_BINS_1D` | 200 | bhattacharyya.py | Bins for 1D marginal BC |
| `N_BINS_2D` | 80 | bhattacharyya.py | Bins per axis for 2D joint BC |
| `SUBSAMPLE` | 10 | kde_contours.py | Use every Nth frame for KDE |
| `GRID_PTS` | 200 | kde_contours.py | KDE evaluation grid resolution |

## Interpreting the Bhattacharyya heatmaps

The Bhattacharyya coefficient (BC) quantifies how much two probability
distributions overlap, ranging from 0 (no overlap) to 1 (identical):

$$
\text{BC}(p, q) = \int \sqrt{\,p(\mathbf{x})\,q(\mathbf{x})\,}\;\mathrm{d}\mathbf{x}
$$

Estimated here via density-normalised histograms.

**Diagonal cells** show the **1D marginal BC** (200 bins) — i.e. how well the
oxidised and reduced projections overlap along each individual PC axis. A value
near 1.0 means the two states sample the same range of that PC; a lower value
means the states are separated along that component.

**Off-diagonal cells** show the **2D joint BC** (80 × 80 bins) for each pair
of PCs. This captures whether the two states occupy the same region of the
joint PC_i–PC_j plane. Even if two PCs overlap individually (high diagonal
values), their 2D joint overlap can be lower if the states differ in how those
PCs are correlated — for example, one state might couple PC1 and PC2
positively whilst the other does not.

In short: the diagonal tells you *which PCs separate the states*, and the
off-diagonal tells you *which pairs of PCs separate them in combination*.

## Reference

Bhattacharyya, A. (1943). On a measure of divergence between two statistical
populations defined by their probability distributions. *Bulletin of the
Calcutta Mathematical Society*, 35, 99–109.
