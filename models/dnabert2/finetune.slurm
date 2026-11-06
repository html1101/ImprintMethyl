#!/bin/bash
# Bridges-2 job for finetune.py. Edit -A and the -o path for your account.
# L40S instead of H100: --gres=gpu:l40s-48:1 and --constraint=l40s (about 2.3x slower).
# Shorter job that starts sooner: -t 0:50:00 and MAXTRAIN=150000 MAXTEST=20000.
#SBATCH -A cis260266p
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:h100-80:1
#SBATCH --constraint=h100
#SBATCH -N 1
#SBATCH -c 8
#SBATCH --mem=64G
#SBATCH -t 2:00:00
#SBATCH -J dnabert2-ft
#SBATCH -o /ocean/projects/cis260266p/rpandit1/imprintmethyl/runs/ft-%j.out
module load anaconda3/2024.10-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate dnabert2
export IMPRINT_BASE=${IMPRINT_BASE:-/ocean/projects/cis260266p/rpandit1/imprintmethyl}
export HF_HOME=$IMPRINT_BASE/hf_cache HF_HUB_OFFLINE=1 TOKENIZERS_PARALLELISM=false
nvidia-smi --query-gpu=name --format=csv,noheader
python $IMPRINT_BASE/dnabert2/finetune.py
