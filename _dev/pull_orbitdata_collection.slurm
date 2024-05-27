#!/bin/bash

#SBATCH --output=../outputs/collection_%a.out
#SBATCH --nodes=1                     
#SBATCH --ntasks=8
#SBATCH --cpus-per-task=2
#SBATCH --time=4:00:00   
#SBATCH --partition=high_priority
#SBATCH --account=gbesla  
#SBATCH --qos=user_qos_gbesla
#SBATCH --array 00-99

### module load python/3.8
source activate tart

date
echo "started pulling orbits for snapshot $SLURM_ARRAY_TASK_ID"
python3.8 pull_orbitdata-collection.py $SLURM_ARRAY_TASK_ID "high" "major"
python3.8 pull_orbitdata-collection.py $SLURM_ARRAY_TASK_ID "high" "minor"
python3.8 pull_orbitdata-collection.py $SLURM_ARRAY_TASK_ID "low" "major"
python3.8 pull_orbitdata-collection.py $SLURM_ARRAY_TASK_ID "low" "minor"

echo "finished"
date


