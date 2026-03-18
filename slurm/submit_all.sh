#!/usr/bin/env bash
# Submit the PCA pipeline with SLURM dependency chaining.
#
# Reads pre-existing Cα-aligned streams from analysis/work/pca/cartesian/.
#
# 02_cart + 03_dihed [parallel]  →  04_figures
#
# Usage:
#   bash PCA/slurm/submit_all.sh

set -euo pipefail
cd /user/work/lh22089

# Clean previous outputs
rm -f PCA/data/*.dat PCA/figures/*.pdf
mkdir -p PCA/data PCA/figures PCA/logs

# Steps 1+2: Cartesian and dihedral PCA (parallel)
JID1=$(sbatch --parsable PCA/slurm/02_cart_pca.slurm)
JID2=$(sbatch --parsable PCA/slurm/03_dihed_pca.slurm)
echo "Submitted 02_cart_pca: $JID1"
echo "Submitted 03_dihed_pca: $JID2"

# Step 3: Python figures (after both PCA jobs complete)
JID3=$(sbatch --parsable --dependency=afterok:${JID1}:${JID2} PCA/slurm/04_figures.slurm)
echo "Submitted 04_figures: $JID3  (after $JID1 + $JID2)"

echo ""
echo "Pipeline: [$JID1 + $JID2] → $JID3"
