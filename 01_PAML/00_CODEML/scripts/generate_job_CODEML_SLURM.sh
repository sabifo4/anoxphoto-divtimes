#!/bin/bash

# Get args
aln=$1         # 1, 2, etc.
pipedir=$2     # Path to pipeline dir
name_wd=$3     # Name of the working directory, e.g., `euk110`
acc_name=$4    # Name of the account to ue
part_name=$5   # Name of the partition to use
job_name=$6    # Name that will be used to tag this job
codeml_exec=$7 # Name of the CODEML binary
bool_paml=$8   # Boolean, TRUE if PAML is in PATH, FALSE otherwise
req_ram=$9     # Requested RAM. E.G., `2G`
nolim_t=$10    # Boolean, TRUE if unlimited time to submit job, FALSE otherwise
# Replace vars in template bash script for job array
cp pipeline_Hessian_CODEML_template_SLURM.sh $pipedir/pipeline_Hessian.sh
# Replace num of job arrays to take place
if [[ $aln -eq 1 ]]
then 
sed -i 's/array\=1..*/array\=1/' $pipedir/pipeline_Hessian.sh
else 
sed -i 's/NUM/'${aln}'/' $pipedir/pipeline_Hessian.sh
fi
# Replace account name and job name
sed -i 's/ACCNAME/'${acc_name}'/' $pipedir/pipeline_Hessian.sh
sed -i 's/JOBNAME/'${job_name}'/' $pipedir/pipeline_Hessian.sh
# Replace partition name
if [[ $part_name =~ "NA" ]]
then
sed -i '/^..*PARTNAME..*$/d' $pipedir/pipeline_Hessian.sh
else
sed -i 's/PARTNAME/'${part_name}'/' $pipedir/pipeline_Hessian.sh
fi
# Replace requested RAM
sed -i 's/REQRAM/'${req_ram}'/' $pipedir/pipeline_Hessian.sh
# Replace name of working directory
upd_wd=$( echo $name_wd | sed 's/\//\\\//g' | sed 's/_/\\_/g' )
sed -i 's/WDNAME/'${upd_wd}'/g' $pipedir/pipeline_Hessian.sh
# Add name for CODEML binary
sed -i 's/CODEMLEXEC/'${codeml_exec}'/g' $pipedir/pipeline_Hessian.sh
# Replace relative path to main dir by just the executable name
if [[ $bool_paml =~ "Y" ]]
then
sed -i 's/\$main_dir\/'${codeml_exec}'/'${codeml_exec}'/' $pipedir/pipeline_Hessian.sh
fi
# If no time allocated, remove line
if [[ ${nolim_t} =~ "Y" ]]
then
sed -i '/^..*time..*Runtime..*$/d' $pipedir/pipeline_Hessian.sh
fi
