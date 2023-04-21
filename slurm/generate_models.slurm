#!/bin/bash
#SBATCH --job-name=generate_models # Job name
#SBATCH --mail-type=END,FAIL         # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=jk4730@columbia.edu    # Where to send mail	
#SBATCH --nodes=1                    # Run all processes on a single node	
#SBATCH --ntasks=1                   # Number of processes
#SBATCH --cpus-per-task=10             # Number of cpus per task
#SBATCH --mem=16gb                    # Total memory limit
#SBATCH --time=06:00:00              # Time limit hrs:min:sec
#SBATCH --output=generate_models_mp_%j.log # Standard output and error log

module purge
eval "$(conda shell.bash hook)"
conda activate myenv

python script.py