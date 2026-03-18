#!/usr/bin/env bash
#SBATCH --job-name=pca_figs
#SBATCH --account=ptch000721
#SBATCH --partition=compute
#SBATCH --time=01:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G
#SBATCH --output=PCA/logs/figures_%j.out
#SBATCH --error=PCA/logs/figures_%j.err

set -euo pipefail
cd /user/work/lh22089
mkdir -p PCA/figures

PYTHON=/user/work/lh22089/miniforge3/bin/python

echo "Bhattacharyya coefficient heatmaps..."
$PYTHON PCA/scripts/python/bhattacharyya.py

echo "KDE contour plots..."
$PYTHON PCA/scripts/python/kde_contours.py

echo "Done:"
ls PCA/figures/
