#!/bin/bash
# Bridges-2 throughput pilot (random 1 kb sequences). Edit -A and the -o path for your account.
#SBATCH -A cis260266p
#SBATCH -p GPU-shared
#SBATCH --gres=gpu:h100-80:1
#SBATCH --constraint=h100
#SBATCH -N 1
#SBATCH -c 4
#SBATCH --mem=32G
#SBATCH -t 20:00
#SBATCH -J dnabert2-pilot
#SBATCH -o /ocean/projects/cis260266p/rpandit1/imprintmethyl/dnabert2/pilot-%j.out
module load anaconda3/2024.10-1
source $(conda info --base)/etc/profile.d/conda.sh
conda activate dnabert2
export HF_HOME=/ocean/projects/cis260266p/rpandit1/imprintmethyl/hf_cache HF_HUB_OFFLINE=1
cd /ocean/projects/cis260266p/rpandit1/imprintmethyl/dnabert2
nvidia-smi --query-gpu=name,memory.total --format=csv
python pilot.py || { echo "=== retry without triton"; NO_TRITON=1 python pilot.py; }
