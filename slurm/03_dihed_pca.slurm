#!/usr/bin/env bash
#SBATCH --job-name=pca_dihed
#SBATCH --account=ptch000721
#SBATCH --partition=compute
#SBATCH --time=12:00:00
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=64G
#SBATCH --output=PCA/logs/dihed_pca_%j.out
#SBATCH --error=PCA/logs/dihed_pca_%j.err

set -euo pipefail
cd /user/work/lh22089
mkdir -p PCA/data

echo "Dihedral PCA: covariance + projections..."
cpptraj -i PCA/scripts/cpptraj/03_dihedral_pca.cpptraj
echo "Done: PCA/data/dihed_*.dat"
