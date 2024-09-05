#!/bin/bash
#SBATCH --chdir=./              # Set to current working directory
#SBATCH --job-name=JOBNAME      # Job name
#SBATCH --account=ACCNAME       # Project code 
#SBATCH --time=15-00:00:00      # Runtime in HH:MM:SS
#SBATCH --partition=PARTNAME    # Partition to use
#SBATCH --mem=REQRAM            # Requested RAM
#SBATCH --array=1-NUM           # Array tasks

#==========================================================================================#
# Contact Sandra Alvarez-Carretero for any doubts about this script: sandra.ac93@gmail.com #
#==========================================================================================#

# ------------------------------------- #
# Creating file structure to run CODEML #
# ------------------------------------- # 

# 1. Find global dirs for paths
pipeline_dir=$( pwd )
main_dir=$( echo $pipeline_dir | sed 's/WDNAME..*/WDNAME/' )
cd $main_dir/Hessian/$SLURM_ARRAY_TASK_ID 
home_dir=$( pwd ) 

# 3. Create specific log file
exec 3>&1> >(while read line; do echo "$line" >> $pipeline_dir/log.hessian.dir$SLURM_ARRAY_TASK_ID .txt; done;) 2>&1
start=`date`
echo Job starts":" $start

# 4. Start analysis
echo The analyses will take place in directory $home_dir
printf "\n"
# Move to analysis dir
cd $home_dir
# Soft link the tmp* files here 
ln -s $home_dir/prepare_codeml/tmp0001.ctl .
ln -s $home_dir/prepare_codeml/tmp0001.trees .
ln -s $home_dir/prepare_codeml/tmp0001.txt .

# 5. Run CODEML
printf "\nRunning CODEML to calculate the Hessian and the gradient...\n"
$main_dir/CODEMLEXEC tmp0001.ctl

# 6. Close
printf "\n"
echo CODEML FINISHED"!"
printf "\n"
end=`date`
echo Job ends":" $end
