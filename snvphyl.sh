#!/bin/bash

#Number of nodes to request. Leave at 1 unless you're running an MPI job.
#SBATCH -N 1
# Set number of processors you want for your job here.
#SBATCH --ntasks=14 
# Set the amount of memory you want (in megabytes)
#SBATCH --mem=40000
# Set amount of time to give your job here. (D-HH:MM)
#SBATCH --time=1-00:00
# Standard out will go here (%j will be replaced by job ID when script is run).
#SBATCH -o /mnt/nas/slurm_logs/job_%j.out
# Standard error will go here (percent j replaced by job ID when script is run).
#SBATCH -e /mnt/nas/slurm_logs/job_%j.err

# Your code goes here. Example activates a virtualenv and then runs a script that needs that virtualenv.
source /mnt/nas/Virtual_Environments/snvphylcli/bin/activate
