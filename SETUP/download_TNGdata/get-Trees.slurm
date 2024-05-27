#!/bin/bash

#SBATCH --output=trees%a.out
#SBATCH --mem=1gb                    
#SBATCH --time=02:00:00   
#SBATCH --partition=standard
#SBATCH --account=gbesla  
#SBATCH --array 0-8
#SBATCH --qos=user_qos_gbesla

cd /xdisk/gbesla/mig2020/rsgrps/gbeslastudents/IllustrisTNG/TNG100-1/trees
wget -nc --no-check-certificate --content-disposition --header="API-Key: <enter API key here>" "http://www.tng-project.org/api/TNG100-1/files/sublink.$SLURM_ARRAY_TASK_ID.hdf5"


