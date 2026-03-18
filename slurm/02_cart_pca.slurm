#!/usr/bin/env bash
#SBATCH --job-name=pca_cart
#SBATCH --account=ptch000721
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=32G
#SBATCH --output=PCA/logs/cart_pca_%j.out
#SBATCH --error=PCA/logs/cart_pca_%j.err

set -euo pipefail
cd /user/work/lh22089
mkdir -p PCA/data

echo "Cartesian PCA: covariance + projections..."
cpptraj -i PCA/scripts/cpptraj/02_cartesian_pca.cpptraj
echo "Done: PCA/data/cart_*.dat"
