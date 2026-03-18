#!/usr/bin/env bash
# Run the PCA pipeline sequentially (non-SLURM fallback).
#
# Reads pre-existing Cα-aligned streams from analysis/work/pca/cartesian/.
#
# Usage:
#   bash PCA/run.sh

set -euo pipefail

PROJECT=/user/work/lh22089
CPPTRAJ=PCA/scripts/cpptraj
PYTHON=$PROJECT/miniforge3/bin/python

cd "$PROJECT"

# Clean previous outputs
rm -f PCA/data/*.dat
rm -f PCA/figures/*.pdf
mkdir -p PCA/data PCA/figures

echo "=== Step 1/4: Cartesian PCA ==="
cpptraj -i $CPPTRAJ/02_cartesian_pca.cpptraj

echo "=== Step 2/4: Dihedral PCA ==="
cpptraj -i $CPPTRAJ/03_dihedral_pca.cpptraj

echo "=== Step 3/4: Bhattacharyya coefficient heatmaps ==="
$PYTHON PCA/scripts/python/bhattacharyya.py

echo "=== Step 4/4: KDE contour plots ==="
$PYTHON PCA/scripts/python/kde_contours.py

echo "=== Done ==="
echo "Data:    PCA/data/"
echo "Figures: PCA/figures/"
ls PCA/figures/
