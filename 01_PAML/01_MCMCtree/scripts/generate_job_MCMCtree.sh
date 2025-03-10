#!/bin/bash

# Get args
dir=$1         # Alignment #1, #2, #3... The alignment being parsed at the moment
clock=$2       # `GBM` or `ILN`
ndat=$3        # 1, 2, 3... As many blocks as partitions in the alignment
pipeloc=$4     # Path to MCMCtree pipeline dir
runmcmc=$5     # Command to execute MCMCtree
nchains=$6     # Number of MCMCs to run
name_wd=$7     # Name working directory
duplication=$8 # Enable duplication option or not

# Replace vars in template bash script for job array
cp pipeline_MCMCtree_template.sh $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/DIR/'${dir}'/g' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/CLK/'${clock}'/g' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/NUMPARTS/'${ndat}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/CMDRUN/'${runmcmc}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/DUP_BOOL/'${duplication}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"

if [[ $nchains -eq 1 ]]
then 
sed -i 's/for\ TASK\_ID\ in..*//' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/^done//' $pipeloc/$dir/$clock/pipeline_$clock".sh"
sed -i 's/^do//' $pipeloc/$dir/$clock/pipeline_$clock".sh"
else 
sed -i 's/NUM/'${nchains}'/' $pipeloc/$dir/$clock/pipeline_$clock".sh"
fi

# Replace name of working directory
upd_wd=$( echo $name_wd | sed 's/\//\\\//g' | sed 's/_/\\_/g' )
sed -i 's/WDNAME/'${upd_wd}'/g' $pipeloc/$dir/$clock/pipeline_$clock".sh"