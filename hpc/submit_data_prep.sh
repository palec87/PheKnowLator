#!/bin/bash
#SBATCH --job-name=pheknowlator_data_prep
#SBATCH --output=logs/data_prep_%j.out
#SBATCH --error=logs/data_prep_%j.err
#SBATCH --account=f202500006hpcvlabualga
#SBATCH --partition=dev-x86
#SBATCH --time=4:00:00
#SBATCH --mem=64G
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4


# PheKnowLator Data Preparation - SLURM Submission Script
# ========================================================
# 
# This script submits the data preparation job to an HPC cluster using SLURM.
# Adjust the resource requirements above based on your cluster's configuration.
#
# Usage:
#   sbatch submit_data_prep.sh [step]
#
# Examples:
#   sbatch submit_data_prep.sh              # Run all steps
#   sbatch submit_data_prep.sh genomic_ids  # Run only genomic identifiers step

# Load required modules (adjust for your HPC environment)
# module load python/3.8
# module load conda

# Activate conda environment if using one
# source activate ontology

# Set up environment variables
export PYTHONUNBUFFERED=1
export OMP_NUM_THREADS=$SLURM_CPUS_PER_TASK

# Create logs directory if it doesn't exist
mkdir -p logs

# activate conda environment
module load Miniconda3
conda activate kgjava 

# Print job information
echo "======================================"
echo "Job ID: $SLURM_JOB_ID"
echo "Job Name: $SLURM_JOB_NAME"
echo "Node: $SLURM_NODELIST"
echo "Start Time: $(date)"
echo "======================================"

# Get step argument if provided
STEP=${1:-all}

# Run the data preparation script
python data_preparation_hpc.py \
    --log-dir ./logs \
    --data-dir ../resources \
    --step $STEP \
    --skip-downloads

# Capture exit code
EXIT_CODE=$?

# Print completion information
echo "======================================"
echo "End Time: $(date)"
echo "Exit Code: $EXIT_CODE"
echo "======================================"

exit $EXIT_CODE
