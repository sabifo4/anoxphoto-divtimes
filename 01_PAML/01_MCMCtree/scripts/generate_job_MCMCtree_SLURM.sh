#!/bin/bash

# Get args
dir=$1            # Alignment #1, #2, #3... The alignment being parsed at the moment
clock=$2          # `GBM` or `ILN`
ndat=$3           # 1, 2, 3... As many blocks as partitions in the alignment
pipeloc=$4        # Path to MCMCtree pipeline dir
runmcmc=$5        # Command to execute MCMCtree
nchains=$6        # Number of MCMCs to run
name_wd=$7        # Name working directory
duplication=$8    # Enable duplication option or not
checkp=$9         # Enable checkpointing?
# Args based on the SLURM job
acc_name=${10}    # Name of the account to ue
part_name=${11}   # Name of the partition to use
job_name=${12}    # Name that will be used to tag this job
bool_paml=${13}   # Boolean, TRUE if PAML is in PATH, FALSE otherwise
req_ram=${14}     # Requested RAM. E.G., `2G`
nolim_t=${15}     # Boolean, TRUE if unlimited time to submit job, FALSE otherwise
tagtime=${16}     # Character, time to run the job

# Replace vars in template bash script for job array
cp pipeline_MCMCtree_template_SLURM.sh $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/DIR/'${dir}'/g' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/CLK/'${clock}'/g' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/NUMPARTS/'${ndat}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/CMDRUN/'${runmcmc}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/DUP_BOOL/'${duplication}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/CHECKP_BOOL/'${checkp}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"

# Add number of chains
# Replace only first encounter
sed -i 's/NUM/'${nchains}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"

# Replace name of working directory
upd_wd=$( echo $name_wd | sed 's/\//\\\//g' | sed 's/_/\\_/g' )
sed -i 's/WDNAME/'${upd_wd}'/g' $pipeloc/$dir/$clock/pipeline_$clock".sh"

# Replace account job name
sed -i 's/JOBNAME/'${job_name}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
# Replace account name
if [[ ${acc_name} =~ "NA" ]]
then
sed -i '/^..*ACCNAME..*$/d' $pipeloc/$dir/$clock/pipeline_$clock".sh"
else
sed -i 's/ACCNAME/'${acc_name}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
fi
# Replace partition name
if [[ ${part_name} =~ "NA" ]]
then
sed -i '/^..*PARTNAME..*$/d' $pipeloc/$dir/$clock/pipeline_$clock".sh"
else
sed -i 's/PARTNAME/'${part_name}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
fi
# Replace requested RAM
sed -i 's/REQRAM/'${req_ram}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
# Replace relative path to main dir by just the executable name
if [[ ${bool_paml} =~ "Y" ]]
then
sed -i 's/\$main_dir\/'${runmcmc}'/'${runmcmc}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
fi
# If time allocated, then fix, otherwise add:
# #SBATCH --time=00:00:00            # Runtime in HH:MM:SS
if [[ ${nolim_t} =~ "Y" ]]
then
sed -i 's/TAGTIME/00\:00\:00/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
elif [[ ${nolim_t} =~ "N" ]]
then
tagtime=$( echo $tagtime | sed 's/\:/\\\:/g' | sed 's/-/\\\-/g' )
sed -i 's/TAGTIME/'${tagtime}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
fi
