# Bayesian inference of species divergences

## 1. Setting the file structure to run `MCMCtree`

Given that we already generated the calibrated tree, alignment, and `in.BV` files... We have everything we need to run `MCMCtree`!

First, we will create the file structure required for the timetree inference analyses using the following code snippet:

```sh
# Run from `spiders_dating` dir on your HPC
# Please change directories until
# you are there
# Then, run the following commands.
num_aln=3     # num_dirs --> 1: "noeurycyde", 2: "noUCEs", 3: "supaln"
num_chains=16 # num chains we will run
for j in `seq 1 $num_chains`
do
for i in `seq 1 $num_aln`
do
mkdir -p MCMCtree/$j/{GBM,ILN}/$i
done
done
# Set pipelines dir
for i in `seq 1 $num_aln`
do
mkdir -p pipelines_MCMCtree/$i/{GBM,ILN}
done
```

The `spiders_dating` directory  will now have these two extra directories with the corresponding subdirectories:

```text
spiders_dating
  |- MCMCtree
  |    |- [1-16] # 16 chains
  |         |- [GBM|ILN]
  |              |- [1-3] # One directory per dataset | 1: "noeurycyde", 2: "noUCEs", 3: "supaln"
  |                
  |- pipelines_MCMCtree
       |- [1-3] # One for each dataset
            |- [GBM|ILN]/
```

Now, we will transfer our in-house bash scripts and corresponding template files with which we will parse our data and generate the job arrays that will be submitted to the HPC to run `MCMCtree`:

```sh
# Run from `01_PAML/01_MCMCtree/scripts`
rsync -avz --copy-links *sh <uname>@<logdetails>:<path>/spiders_dating/scripts
```

You can now go back to your HPC and run the [`generate_job_MCMCtree.sh` script](scripts/generate_job_MCMCtree.sh) (script that should have already been transferred) by using the commands below:

```sh
# Run from the `spiders_dating` dir on your HPC
# Please change directories until you are there
# Then run the following commands
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_aln=3 # num alignment
num_chains=16
# Arg 1: Number of gene alignments being analysed
# Arg 2: Clock model (e.g., "GBM", "ILN", or "CLK)
# Arg 3: Number of partitions in the alignment file
# Arg 4: Path to the pipeline directory
# Arg 5: Command to execute MCMCtree (e.g., `mcmctree`, `mcmctree_2000`, etc.)
# Arg 6: Number of chains run
# Arg 7: Name of working directory (e.g., `spiders_dating`)
# Arg 8: Boolean, enable duplication option? 0: no, 1: yes
# 
for i in `seq 1 $num_aln`
do
./generate_job_MCMCtree.sh $i GBM 1 $home_dir/pipelines_MCMCtree mcmctree4.10.7 $num_chains spiders_dating 0
./generate_job_MCMCtree.sh $i ILN 1 $home_dir/pipelines_MCMCtree mcmctree4.10.7 $num_chains spiders_dating 0
done
# Check paths were correctly incorporated!
cd ../pipelines_MCMCtree
grep '^dir=\$' */*/*sh
grep 'MCMCtree' */*/*sh
grep 'alignments' */*/*sh
grep 'Hessian' */*/*sh
```

Now, before running `MCMCtree` to sample from the posterior, we will run `MCMCtree` but sampling from the prior. When doing so, we can check whether the calibration densities (also known as "user-specified priors") constraining some node ages in our fixed tree topology are in conflict with the corresponding marginal densities (i.e., the so-called "effective priors", which the dating program will use to generate the joint time prior together with additional priors such as the birth-death process with sampling and the calibration densities). Oftentimes, truncation issues may arise when the effective priors are in disagreement with the user-specified priors (see an extensive study about this effect in [dos Reis et al. 2015](https://pubmed.ncbi.nlm.nih.gov/26603774/)). To assess the effect of truncation prior to analysing our dataset, we will run `MCMCtree` without the data (i.e., `MCMCtree` will be sampling from the prior distribution, and thus sequence data will not be used).

First, we will generate a directory where `MCMCtree` will run when sampling from the prior:

```sh
# Run from `spiders_dating` dir on your HPC
# Please change directories until you are there
# Then run the following commands
num_aln=3
num_chains=6
for j in `seq 1 $num_chains`
do
for i in `seq 1 $num_aln`
do
mkdir -p MCMCtree_prior/$j/CLK/$i
done
done
```

> [!IMPORTANT]
> When sampling from the prior, the likelihood is not being calculated or estimated. In other words, the most time-consuming part of the MCMC does not take place. To this end, you should be able to gather enough samples with fewer runs than those needed when sampling from the posterior.

Then, we will copy the directory `pipelines_MCMCtree` and will generate a copy called `pipelines_MCMCtree_prior` :

```sh
# Run from `spiders_dating` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
cp -R pipelines_MCMCtree pipelines_MCMCtree_prior
cd pipelines_MCMCtree_prior
```

We will modify the bash script that will be submitted as a job array so that the `userdata` option in the control file is equal to `0`, which enables `MCMCtree` to sample from the prior instead of sampling from the posterior (i.e., the alignment file is ignored). The rest of the setting concerning the evolutionary model will not be enabled by `MCMCtree` as the sequence data are not being used. Lastly, we will change the path to where the results will be stored:

```sh
# Run from `pipelines_MCMCtree_prior` dir on your HPC
# Please change directories until you are there
# Then run the following commands
# 
# Prepare directories to sample from the prior,
# only one needed as nucleotide subsitution models 
# are not used.
#
num_aln=3
home_dir=$( pwd )
for i in `seq 1 $num_aln`
do
cd $home_dir/$i
rm -r ILN 
mv GBM CLK
cd CLK
mv pipeline_GBM.sh pipeline_CLK.sh
done
cd $home_dir

# Modify bash script: options `usedata` and `model`
sed -i "s/..*usedata..*/sed \-i \'s\/usedata\.\.\*\/usedata \= 0\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"\nsed \-i \'s\/model\.\.\*\/model \= 4\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"/" */*/*sh
# Modify clock model 
sed -i "s/^fi/fi\nsed \-i \'s\/clock\.\.\*\/clock \= 1\/\' \$home\_dir\/mcmctree\_\$dir\"_r\"\$SGE\_TASK\_ID\"\.ctl\"/" */*/*sh

# Modify path to save results 
sed -i 's/\/GBM\//\/CLK\//' */*/*sh
sed -i 's/MCMCtree/MCMCtree\_prior/' */*/*sh 
# Comment soft link
sed -i 's/ln \-s/\#ln \-s/' */*/*sh
# Change number of chains to 6!
sed -i 's/-t 1-16/-t 1-6/' */*/*sh
```

Now, you can check that the lines have been correctly modified:

```sh
# Run from `pipelines_MCMCtree_prior` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
cd ../
grep '^dir=' */*/*sh
grep 'usedata' */*/*sh
grep 'model'  */*/*sh
grep 'MCMCtree' */*/*sh
grep '#$ -t' */*/*sh
```

## 2. Analyses with `MCMCtree` when sampling from the prior

### Run `MCMCtree` in the HPC - prior

Now, we will be able to run `MCMCtree` first when sampling from the prior (i.e., no data used!) using the code snippet below:

```sh
# Run from `pipelines_MCMCtree_prior/CLK` dir on your HPC.
# Please change directories until
# you are there. Then run the following
# commands.
chmod 775 *sh
qsub pipeline_CLK.sh
```

### Setting the file structure to analyse `MCMCtree` output - prior

We will now create a `sum_analyses` directory to analyse the `MCMCtree` output. Nevertheless, we first need to transfer the data from the cluster to the corresponding directory on our local PC for further analyses:

```sh
# Run everything from `spiders_dating` in your HPC
num_chains=6
num_datasets=3 # 1: "noeurycyde", 2: "noUCEs", 3: "supaln"
mkdir -p tmp_to_transfer/00_prior
cd tmp_to_transfer
for j in `seq 1 $num_datasets`
do
for i in `seq 1 $num_chains`
do
mkdir -p 00_prior/CLK/$j/$i/
# Now, copy the files that are required for sum stats
# We have run 6 chains for analyses sampling from the prior
printf "\n[[ Copying run "$i" for analyses sampling from the prior -- dataset "$j" ]]\n\n"
cp ../MCMCtree_prior/$i/CLK/$j/mcmc.txt 00_prior/CLK/$j/$i
cp ../MCMCtree_prior/$i/CLK/$j/*ctl 00_prior/CLK/$j/$i
cp ../MCMCtree_prior/$i/CLK/$j/SeedUsed 00_prior/CLK/$j/$i
done
grep 'Species tree for FigTree' -A1 ../MCMCtree_prior/$j/CLK/1/out.txt | sed -n '2,2p' > 00_prior/node_tree_$j.tree
done
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `01_PAML/01_MCMCtree` dir on your local PC
# Please change directories until you are there
# Then run the following commands
# If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
mkdir sum_analyses
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/spiders_dating/tmp_to_transfer/00_prior .
rsync -avz --copy-link <uname>@<logdetails>:<path>/spiders_dating/pipelines_MCMCtree_prior .
# Remove blank output files
rm pipelines_MCMCtree_prior/*/*/*sh.o*
```

### MCMC diagnostics - prior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to run MCMC diagnostics!

First, however, we need to generate a file with calibration information that is compatible with the subsequent scripts. For that purpose, we can use our in-house R script [`Merge_node_labels.R](scripts/Merge_node_labels.R), which will generate one calibration file for each dataset analysed in case there are differences in the tree topologies and the node labels which age is being constrained.

---

> [!NOTE]
> In order to incorporate the node names in the summary files that will be generated in the subsequent steps, we have created a copy of the calibration files that the script [`Merge_node_labels.R](scripts/Merge_node_labels.R) will have generated (i.e., comma-separated files that you shall find now saved in directory [calib_files](calib_files)) and added additional node names using the code snippet below:

```sh
# Run from `01_MCMCtree`
cd calib_files
for i in *csv
do
name=$( echo $i | sed 's/\.csv/\_margVScalib\.csv/' )
cp $i $name
if [[ $name =~ "noeurycyde" ]]
then
printf "Ascorhynchidae;334;0\n" >> $i
fi
printf "Nymphonoidea;240;0\n" >> $i
printf "Ausctrodecidae;290;0\n" >> $i
printf "Pycnogonidae;295;0\n" >> $i
printf "Endeidae;350;0\n" >> $i
printf "Phoxichilidiidae;357;0\n" >> $i
printf "Pallenopsidae;374;0\n" >> $i
printf "Ammotheidae;392;0\n" >> $i
done
```

---

Now, we can run the R script [`MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R) and follow the detailed step-by-step instructions detailed in the script. In a nutshell, the protocol will be the following:

1. Load the `mcmc.txt` files generated after each run.
2. Generate a convergence plot with the unfiltered chains.
3. Find whether there are major differences between the time estimates sampled across the chains for the same node sin the 97.5% and the 2.5% quantiles. If so, flag and delete said chains.
4. If some chains have not passed the filters mentioned above, create an object with the chains that have passed the filters.
5. Generate a new convergence plot with those chains that passed the filters.
6. Calculate Rhat, tail-ESS, and bulk-ESS to check whether chain convergence has been reached with the chains that have passed filters.

The MCMC diagnostics did not find any of the chains problematic after running [our in-house R script `MCMC_diagnostics_prior.R`](scripts/MCMC_diagnostics_prior.R). Therefore, we used [our in-house bash script `Combine_MCMC.sh`](scripts/Combine_MCMC.sh) to concatenate all the `mcmc.txt` files for the 6 chains in a unique file.

```sh
# Run from `01_MCMCtree/scripts`
cp Combine_MCMC.sh ../sum_analyses/00_prior
# One argument taken: number of chains
cd ../sum_analyses/00_prior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too
dataset=$( echo CLK )
num_dat=3
name_dat=( 'noeurycyde' 'noUCEs' 'supaln' )
count=-1 #  counter will start at 0 in first iteration!
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
./Combine_MCMC.sh $dataset/$i mcmc_files_${name_dat[count]}_CLK "`seq 1 6`" CLK 20000 Y ${name_dat[count]}_CLK
done
```

The script above will generate directories called `mcmc_files*_CLK` inside the `00_prior` directory, where the `mcmc.txt` with the concatenated samples will be saved. In addition, directories with individual `mcmc.txt` files of those chains that passed the filters will be created (i.e., see `mcmcf4traces*_CLK` directories); you can read such files in programs like `Tracer` to assess the traces and run other visual MCMC diagnostics.

We will now create a dummy alignment with only 2 AAs, which is required to generate the `FigTree` files with the mean time estimates obtained when using the concatenated `mcmc.txt` files. In order to do that, we can run the `Generate_dummy_aln.R`. Once you run this script, a new directory called `dummy_aln` will be created, which will contain the input dummy alignment.

We have also generated dummy control file to read the dummy alignment. Additionally, we have enabled option `print = -1`. This print setting lets `MCMCtree` know that an MCMC is not to be run. Instead, `MCMCtree` is told to read the input files (file with the dummy alignment, the calibrated tree file, and the concatenated `mcmc.txt` file) and summarise the samples in the `mcmc.txt` (those that were collected from those chains that passed the filters!). The final mean estimated divergence times and the corresponding CIs will be written in the output `FigTree.tre` file.

```sh
##> [IMPORTANT] Before running the `for` loop below,
##> please check the commented sections and ammend the
##> commands accordingly depending on whether you are
##> using a program that is exported to the system's path
##> or a binary that needs to be execute with a relative
##> path

# Run from `sum_analyses/00_prior`
name_dat=( 'noeurycyde' 'noUCEs' 'supaln' )
num_dat=3
count=-1
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
printf "\n[[ Analysing dataset "${name_dat[count]}" ]]\n"
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls ${name_dat[count]}"_calib_MCMCtree.tree"`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_${name_dat[count]}"_CLK"
printf "[[ Generating tree file for concatenated \"mcmc.txt\"  ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_CLK_95HPD.tree"
printf "\n"
##> [IMPORTANT] If you have PAML v4.9h, you can also obtain the timetree with 95%CIs
##> Modify the command below if you have a different alias to run `MCMCtree` from
##> PAML v4.9h (here we are using `mcmctree49h_sum95CI`) than
##> the one used below and uncomment the commands below
##> If not, please just ignore and run the `for` loop as is!
#
#tmp_tt=$( echo tree_${name_dat[count]}"_uncalib.tree" )
#cp $tt_dir/$tmp_tt dummy_cal.tree
#cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
#sed -i 's/treefile..*/treefile\ \=\ dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
#sed -i "s/\;/\'B\(5\.14\,6\.361\,1e\-300\,0\.025\)\'\;/" dummy_cal.tree
#mcmctree49h_sum95CI *ctl mcmctree_dummy_95CI.ctl
#printf "\n"
#mv FigTree.tre FigTree_${name_dat[count]}"_CLK_95CI.tree"
cd $base_dir
done
```

The next step is to plot the user-specified prior VS the effective prior. We used our in-house R script [`Check_priors_effVSuser.R`](scripts/Check_priors_effVSuser.R) to generate these plots. If you are to run this script with other datasets, however, make sure that your "hard bounds" are not `0.000` in the `Calibnodes_*csv` files and, instead, they are `1e-300` (i.e., while 1e-300 is rounded to `0.000` in the `MCMCtre` output, which can be used to generate the csv files aforementioned, we need `1e-300` to plot distributions in R). To make sure this was not affecting our csv files, we ran the following code snippet:

```sh
# Run from `01_MCMCtree/calib_files`
sed -i 's/0\.000/1e\-300/g' *csv
```

The next step is to plot the marginal densities VS the calibration densities. We used our in-house R script [`Check_priors_margVScalib.R`](scripts/Check_priors_margVScalib.R) to generate these plots. Once this script has finished, you will see that a new directory `plots/margVScalib` will have been created. Inside this directory, you will find one directory for each individual dataset with individual plots for each node. In addition, all these plots have been merged into a unique document as well (note: some plots may be too small to see for each node, hence why we have generated individual plots).

Now, once the MCMC diagnostics have finished, you can extract the relevant output that we used to write our manuscript:

```sh
# Run from `01_MCMCtree`
mkdir sum_files_prior
cp -R sum_analyses/00_prior/mcmc_files*CLK/*CLK*tree sum_files_prior/
cp -R sum_analyses/00_prior/CLK/*/*/*all_mean*tsv sum_files_prior/
cp -R plots/ESS_and_chains_convergence/*prior*pdf sum_files_prior/
cp -R plots/margVScalib sum_files_prior/
```

## 3. Analyses with `MCMCtree` when sampling from the posterior

### Run `MCMCtree` in the HPC - posterior

Now that we have verified that there are no issues between the calibration and marginal densities, we can run `MCMCtree` when sampling from the posterior. We will do these analyses under the GBM and ILN relaxed-clock models using the code snippet below:

```sh
# Now, go to directory `spiders_dating/pipelines_MCMCtree/GBM` dir on your HPC
# and run the following command
# Please change directories until you are there.
chmod 775 *sh
qsub pipeline_GBM.sh

# Now, go to directory `pipelines_MCMCtree/ILN` 
# and run the following command
# Please change directories until you are there
chmod 775 *sh
qsub pipeline_ILN.sh
```

> [!IMPORTANT]
> When sampling from the posterior, the likelihood is being calculated or approximated, depending on the `userdata` option you set in the control file to run `MCMCtree`. In other words, the larger the dataset, the more time it will take for `MCMCtree` to finish.

### Setting the file structure to analyse `MCMCtree` output - posterior

We will now create a directory inside the `sum_analyses` directory to analyse the `MCMCtree` output. Nevertheless, we first need to transfer the data from the cluster to the corresponding directory on our local PC for further analyses:

```sh
# Go to your HPC and copy the files that are required for sum stats
# We have run 16 chains for analyses sampling from the posterior
# Therefore, `i` will go form 1 to 16
# Run from `spiders_dating`
cd tmp_to_transfer
num_chains=16
num_datasets=3 # 1: "noeurycyde", 2: "noUCEs", 3: "supaln"
# The `01_posterior` should already exist from the previous analyses
# If not, it will be created during the `for` loop
for j in `seq 1 $num_datasets`
do
for i in `seq 1 $num_chains` 
do
mkdir -p 01_posterior/{GBM,ILN}/$j/$i
printf "\n[[ Copying run "$i" for analyses sampling from the posterior -- dataset "$j" ]]\n\n"
cp ../MCMCtree/$i/GBM/$j/mcmc.txt 01_posterior/GBM/$j/$i
cp ../MCMCtree/$i/GBM/$j/*ctl* 01_posterior/GBM/$j/$i
cp ../MCMCtree/$i/GBM/$j/SeedUsed 01_posterior/GBM/$j/$i
cp ../MCMCtree/$i/ILN/$j/mcmc.txt 01_posterior/ILN/$j/$i
cp ../MCMCtree/$i/ILN/$j/*ctl 01_posterior/ILN/$j/$i
cp ../MCMCtree/$i/ILN/$j/SeedUsed 01_posterior/ILN/$j/$i
done
done
```

Now, you can transfer the temporary directory to the local PC, e.g., using `rsync`:

```sh
# Run from `01_PAML/01_MCMCtree` dir on your local PC
# Please change directories untilyou are there
# Then run the following commands
# If you are running this code with your
# own analyses, make sure that you have correctly
# defined `num_aln` and `num_chains` variables with
# the correct values!
# Note that we will generate some directories for
# when the analyses when sampling from the posterior
# are ready!
cd sum_analyses
# Now, trasnfer the data from the HPC
rsync -avz --copy-links <uname>@<logdetails>:<path>/spiders_dating/tmp_to_transfer/01_posterior .
rsync -avz --copy-links <uname>@<logdetails>:<path>/spiders_dating/pipelines_MCMCtree .
# Remove blank output files
rm pipelines_MCMCtree/*/*/*sh.o*
```

### MCMC diagnostics - posterior

Now that we have the output files from the different MCMC runs in an organised file structure, we are ready to check the chains for convergence!

We are going to run the R script [`MCMC_diagnostics_posterior.R`](scripts/MCMC_diagnostics_posterior.R) and follow the detailed step-by-step instructions detailed in the script, which are essentially the same ones we used when analysing the chains when sampling from the prior. Given that no problems have been found with any of the chains we ran, we are ready to concatenate the parameter values sampled across the 16 independent chains we ran:

```sh
# Run from `01_MCMCtree/scripts`
cp Combine_MCMC.sh ../sum_analyses/01_posterior
# One argument taken: number of chains
cd ../sum_analyses/01_posterior
## Variables needed
## arg1 --> name of directory where analyses have taken place (e.g., CLK, GBM, ILN)
## arg2 --> output dir: mcmc_files_CLK, mcmc_files_GBM, mcmc_files_ILN, etc.
## arg3 --> "`seq 1 36`", "1 2 5", etc. | depends on whether some chains were filtered out or not
## arg4 --> clock model used: ILN, GBM, CLK
## arg5 --> number of samples specified to collect in the control file to run `MCMCtree`
## arg6 --> 'Y' to generate a directory with files compatible with programs such as `Tracer` to visually
##          inspect traceplots, convergence plots, etc. 'N' otherwise
## arg7 --> if arg6 is 'Y', arg7 needs to have a name for the `mcmcf4traces_<name>` that will be
##          generated. If `arg6` is equal to 'N', then please write `N` too
dirname_1=GBM
dirname_2=ILN
num_dat=3
name_dat=( 'noeurycyde' 'noUCEs' 'supaln' )
count=-1 #  counter will start at 0 in first iteration!
for i in `seq 1 $num_dat`
do
count=$(( count + 1 ))
./Combine_MCMC.sh $dirname_1/$i mcmc_files_${name_dat[count]}_GBM "`seq 1 16`" GBM 20000 Y ${name_dat[count]}_GBM
./Combine_MCMC.sh $dirname_2/$i mcmc_files_${name_dat[count]}_ILN "`seq 1 16`" ILN 20000 Y ${name_dat[count]}_ILN
done
```

Once the scripts above have finished, a new directory called `mcmc_files_part_[GBM|ILN]` will be created inside `01_posterior/`, respectively. To map the mean time estimates with the filtered chains, we need to copy a control file, the calibrated Newick tree, and the dummy alignment we previously generated when analysing the results when sampling from the prior:

```sh
##> [IMPORTANT] Before running the `for` loop below,
##> please check the commented sections and ammend the
##> commands accordingly depending on whether you are
##> using a program that is exported to the system's path
##> or a binary that needs to be execute with a relative
##> path

# Run from `sum_analyses_prot/01_posterior` directory.
# Please change directories until
# you are there. Then run the following
# commands.
name_dat=( 'noeurycyde' 'noUCEs' 'supaln' )
num_dat=3
count=-1
for i in `seq 1 $num_dat`
do
printf "\n[[ ANALYSING DATASET "${name_dat[count]}" ]]\n\n"
count=$(( count + 1 ))
base_dir=$( pwd )
cd ../../dummy_ctl_files
ctl_dir=$( pwd )
cd ../../../00_data_formatting/01_inp_data/
tt_dir=$( pwd )
name_tt=`ls ${name_dat[count]}"_calib_MCMCtree.tree"`
cd $ctl_dir
cd ../dummy_aln
aln_dir=$( pwd )
name_aln=`ls *aln`
cd $base_dir
cd mcmc_files_${name_dat[count]}"_GBM"
printf "[[ Generating tree file for concatenated \"mcmc.txt\" for GBM ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_GBM_95HPD.tree"
tmp_tt=$( echo tree_${name_dat[count]}"_uncalib.tree" )
cp $tt_dir/$tmp_tt dummy_cal.tree
cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
sed -i 's/treefile..*/treefile\ \=\ dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
sed -i "s/\;/\'B\(5\.14\,6\.361\,1e\-300\,0\.025\)\'\;/" dummy_cal.tree
mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
mv FigTree.tre FigTree_${name_dat[count]}"_GBM_95CI.tree"
printf "\n"
cd $base_dir/mcmc_files_${name_dat[count]}"_ILN"
printf "[[ Generating tree file for concatenated \"mcmc.txt\" in "$data"/"$i" for ILN ... ... ]]\n"
cp $ctl_dir/*ctl .
name_mcmc=`ls *mcmc.txt`
sed_aln=$( echo $aln_dir"/"$name_aln | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed_tt=$( echo $tt_dir"/"$name_tt | sed 's/\//\\\//g' | sed 's/_/\\_/g' |  sed 's/\./\\\./g' )
sed -i 's/MCMC/'${name_mcmc}'/' *ctl
sed -i -e 's/ALN/'${sed_aln}'/' *ctl
sed -i 's/TREE/'${sed_tt}'/' *ctl
# Run now MCMCtree after having modified the global vars according to the path to these files
# Then, rename the output tree file so we can easily identify later which tree belongs to
# which dataset easier
mcmctree *ctl
mv FigTree.tre FigTree_${name_dat[count]}"_ILN_95HPD.tree"
printf "\n"
##> [IMPORTANT] If you have PAML v4.9h, you can also obtain the timetree with 95%CIs
##> Modify the command below if you have a different alias to run `MCMCtree` than
##> the one used below and uncomment the commands below
##> If not, please just ignore and run the `for` loop as is!
#
#tmp_tt=$( echo tree_${name_dat[count]}"_uncalib.tree" )
#cp $tt_dir/$tmp_tt dummy_cal.tree
#cp mcmctree_dummy.ctl mcmctree_dummy_95CI.ctl
#sed -i 's/treefile..*/treefile\ \=\ dummy\_cal\.tree/' mcmctree_dummy_95CI.ctl
#sed -i "s/\;/\'B\(5\.14\,6\.361\,1e\-300\,0\.025\)\'\;/" dummy_cal.tree
#mcmctree49h_sum95CI mcmctree_dummy_95CI.ctl
#mv FigTree.tre FigTree_${name_dat[count]}"_ILN_95CI.tree"
#printf "\n"
cd $base_dir
done
```

Now, once the MCMC diagnostics have finished, we can run our [in-house R script](scripts/Check_priors_VS_posteriors.R) to plot the posterior distributions against the prior distributions, which can help better assess how informative the data are and whether there are any serious contradictions between the prior and the posterior distributions.

Lastly, you can extract the relevant output that we used to write our manuscript as it follows:

```sh
# Run from `01_MCMCtree`
mkdir sum_files_post
cp -R sum_analyses/01_posterior/mcmc_files_*/FigTree*tree sum_files_post/
cp -R sum_analyses/01_posterior/*/*/*/*all_mean*tsv sum_files_post/
cp -R plots/priorVSpost*pdf sum_files_post/
cp -R plots/ESS_and_chains_convergence/*post*pdf sum_files_post/
```
