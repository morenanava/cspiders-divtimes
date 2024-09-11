# `BASEML` analysis

## 1. Pick rate prior

We will use a vague gamma distribution centered on a mean evolutionary rate estimated by considering the tree height (molecular distance in substitutions per site) and the mean age for the root of our phylogenies (time unit = 100 Mya). As [our rooted phylogenies](../../00_data_formatting/00_raw_data/trees) have information about the branch lengths, we can use [our R in-house script](scripts/calculate_rateprior.R) to calculate the tree heights. We also have a calibration to constrain the root age: a minimum age of 514 Ma and a maximum age of 636.1 Ma, which average we will use as an approximate age for the root of our phylogenies to estimate the mean evolutionary rate.

By setting a vague shape ($\alpha=2$) for the gamma distribution that we will use as a rate prior, we can account for the uncertainty surrounding our mean rate estimate. If we had more knowledge on the mean rate, however, we should use a narrower prior with a larger $\alpha$ that better represents our prior information.

Now, we have all the information we need to calculate the $\beta$ parameters for the gamma distributions. We have written the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R) to carry out all the tasks mentioned above. You can open this file in RStudio to find out the value of $\beta$ and plot the final prior on the rates. A summary of what you will find in the script is described below:

```text
First, we know that the molecular distance (tree height, distance from the root to present time) is equal to the mean evolutionary rate (in substitutions per site per time unit) times the age of the divergence time at the root (in time unit, which we can define later). If we have estimated our phylogeny, and therefore have estimated the value of the branch lengths, we will be able to calculate the tree height. The units of the tree height will be the following:

tree_height = rate * root_age --> units_tree_height = subst/site/time * time = subst/site

There are various functions we can use to calculate the tree heigt. We have chosen the R function `phytools::nodeHeights`. The maximum height calculated by this function corresponds to the length from the root to the heighest tip.

After calculating the tree height of our phylogeny (in subst/site) and considering the age of the root based on the fossil record or geological events (time unit = 1 Ga = 100 Mya = 1e9 years), we can get a rough estimate for the mean rate:

Time unit = 100 Mya --> mean_rate = tree_height / root_age = (subst/site) / (Mya) = subst/site/Mya (time unit = 100 Mya) 

We also know that the mean of the gamma distribution that we want to use as rate prior is our parameter of interest: the mean evolutionary rate. Therefore:

mean_Gamma = mean_rate = alpha / beta 
Time unit = 100 Mya: mean_rate = alpha / beta --> beta = alpha / mean_rate = 2 / mean_rate

The calibrated tree needs to incorporate the age constraints in the same time unit used to infer the mean evolutionary rate and establish the rate prior (i.e., do not forget to scale the calibrations accordingly if needed!). 
```

If you run the [R script `calculate_rateprior.R`](scripts/calculate_rateprior.R), you will see how all the steps described above take place and a new PDF file with the prior distribution to be used will be generated in a new directory called `out_RData`.

We have then updated our [template control files](control_files) with the $\alpha$ and $\beta$ parameters (as defined using the R script above) for the gamma distributions as rate priors. Note that several options in this control file will be subsequently modified to fit the analysis with this dataset (i.e., you will see some options that have flags in capital letters, which will be replaced with the correct value for said option). Given how shallow this tree is, the clock may be expected to be seriously violated, and thus we have fixed a mean for the `sigma2` parameter (i.e., variation in the clock) as 0.1 using a gamma prior with $\alpha=1$ and $\beta=10$: `sigma2_gamma 1 10​`.

## 2. Set up the file structure

Before running `MCMCtree` using the approximate likelihood calculation to speed up timetree inference, we first need to calculate the vector of branch lengths, the gradient (vector), and the Hessian (matrix). We will use `BASEML` for this purpose as our dataset is a nucleotide alignment.

The file structure we will use is the following:

```text
spiders_dating/
  |
  |- alignments/
  |    |- X/ # Directory for alignment X -- have as many directories as alignments
  |       
  |- control_files/ # Pre-defined control file with flags to be later replaced with specific settings
  |
  |- Hessian/
  |    |- X # Directory for alignment X -- have as many directories as alignments
  |          
  |- pipelines_Hessian # Directory where the pipeline to run `BASEML` will be executed
  |
  |- scripts # Scripts used to prepare control files to run `BASEML
  |
  |- trees
      |- calibrated   # Directory with the calibrated tree for `MCMCtree`
      |- uncalibrated # Directory with the uncalibrated tree for `BASEML`
```

To create the `spiders_dating` file structure, we run the following commands from the PC before transferring to the HPC:

```sh
# Run the following commands from 
# directory `00_BASEML`
mkdir -p HPC/spiders_dating
cd HPC/spiders_dating 
num_aln=3
for i in `seq 1 $num_aln`
do
mkdir -p alignments/$i
mkdir -p Hessian/$i/prepare_baseml
mkdir -p trees/{uncalibrated/$i,calibrated/$i}
mkdir -p control_files/$i
done
mkdir -p pipelines_Hessian
mkdir scripts
```

Once the file structure is created, we can now populate it with the input files we have generated some minutes ago: alignment files, tree files, and control files. You can transfer the files to the HPC as you prefer (below, we show an example of how to do this with `rsync`, a system that we will keep as an example throughout the tutorial):

```sh
# Run from `HPC/spiders_dating`
# Copy alignment
cp ../../../../00_data_formatting/01_inp_data/aln_noeurycyde.phy alignments/1
cp ../../../../00_data_formatting/01_inp_data/aln_noUCEs.phy alignments/2
cp ../../../../00_data_formatting/01_inp_data/aln_supaln.phy alignments/3
# Now, transfer calibrated trees
cp ../../../../00_data_formatting/01_inp_data/noeurycyde_calib_MCMCtree.tree trees/calibrated/1
cp ../../../../00_data_formatting/01_inp_data/noUCEs_calib_MCMCtree.tree trees/calibrated/2
cp ../../../../00_data_formatting/01_inp_data/supaln_calib_MCMCtree.tree trees/calibrated/3
# Transfer uncalibrated trees
cp ../../../../00_data_formatting/01_inp_data/tree_noeurycyde_uncalib.tree trees/uncalibrated/1
cp ../../../../00_data_formatting/01_inp_data/tree_noUCEs_uncalib.tree trees/uncalibrated/2
cp ../../../../00_data_formatting/01_inp_data/tree_supaln_uncalib.tree trees/uncalibrated/3
# Next, copy control files
cp ../../control_files/prepbaseml_noeurycyde.ctl control_files/1
cp ../../control_files/prepbaseml_noUCEs.ctl control_files/2
cp ../../control_files/prepbaseml_supaln.ctl control_files/3
# Last, copy the in-house bash scripts with our pipeline
cp ../../scripts/*sh scripts/
# Once everything is ready, you can transfer this directory to your cluster!
# One way of doing this is by using `rsync`, but you may use other approaches.
# Below, you will find an example of the `rsync` commands you should run once
# you replace the tags with your own credentials.
# First, move one dir back so you are inside `HPC`
cd ../
rsync -avz --copy-links spiders_dating <uname>@<server>:<path_to_your_wd_in_HPC>
```

----

**IMPORTANT INFORMATION REGARDING THE PAML VERSION INSTALLED ON OUR HPC SERVER**
We have compiled the PAML programs available for the latest version of this software in our HPC server (v4.10.7) given that the cross-bracing approach is implemented. There are two ways that you could follow to use PAML software on your HPC:

* If you have an older version exported to your path variable and do not want to change it: note that you will need to compile `MCMCtree` after modifying the source code: `NS` needs to be increased as there are more than 500 taxa in our alignment! In that way, please open the `mcmctree.c` file inside `src` directory and change the line where `NS` is defined to the following: `#define NS            2000`. Then, please save the changes and compile `MCMCtree` following the PAML installation guidelines given on [the PAML GitHub repository](https://github.com/abacus-gene/paml/wiki#installation). Once the software is compiled, please rename the `mcmctree` and `baseml` binaries to `mcmctree4.10.7` and `baseml4.10.7`, respectively, and save them in your `spiders_dating` working directory so that they are launched using relative paths.
* If you want to export the latest PAML version to your path variable: note that you will need to compile `MCMCtree` after modifying the source code: `NS` needs to be increased as there are more than 500 taxa in our alignment! In that way, please open the `mcmctree.c` file inside `src` directory and change the line where `NS` is defined to the following: `#define NS            2000`. Then, please export the path to the `bin` directory where you have all the PAML binaries. You will be able to launch these programs by typing `mcmctree` and `baseml` from any location in your file structure if the path has been properly exported.

Please note that `MCMCtree` and `BASEML` are the PAML programs that we will use during all the steps of timetree inference. All inference analyses are therefore run by calling these programs via relative paths (if you chose the former) or absolute paths (if you chose the latter). Given that `conda` does not yet have the latest PAML version available, these are the two possibilities you have to work with the latest PAML version on your HPC.

----

Now, we need to generate other input files to estimate the Hessian and the gradient: the input control files for `BASEML`. To do this in a reproducible manner, you can use the [script `generate_prepbaseml.sh`](scripts/generate_prepbaseml.sh), which you can find in the [`01_PAML/00_BASEML/scripts`](01_PAML/00_BASEML/scripts) and which you should have just transferred to your HPC. Now, connect to your server and run the next code snippet, where you will execute this script. Specifically, the [`generate_prepbaseml.sh` script](scripts/generate_prepbaseml.sh) needs one argument: the directory name where the alignments are saved: `1`, `2`, and `3` in our case!

```sh
# Run from `spiders_dating/scripts` in the HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
chmod 775 *sh
# In this case, there are three alignments, so
# we can execute our script within a loop
num_aln=3
for i in `seq 1 $num_aln`
do
./generate_prepbaseml.sh $i
done
```

To make sure that all the paths have been properly extracted, you can run the following code snippet:

```sh
# Run from `spiders_dating/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */prepare_baseml/*ctl
grep 'treefile' */prepare_baseml/*ctl
```

## 3. Run `BASEML`

### Preparing input files

Now that we have the input files (alignment and tree files) and the instructions to run `BASEML` (control file) in our HPC server, we will be manually running `MCMCtree` inside each `prepare_baseml` directory (see file structure above) in a special mode that launches `BASEML` for the sole purpose we want: to infer the vectors and matrix required to approximate the likelihood calculation.

```sh
# Run `MCMCtree` from
# `spiders_dating/Hessian/1/prepare_baseml`
# dir on the HPC. 
# Please change directories until
# you are in there.
# The first command to change directories 
# will work if you are still in 
# `main/Hessian`, otherwise ignore and 
# move to such directory with the command
# that best suits your current directory.
# Keep changing the value of "dir" to 1, 2, 3:
dir=1
cd $dir/prepare_baseml
../../../mcmctree4.10.7 prepbaseml*ctl
```

First, you will see that `MCMCtree` starts parsing the first locus. Then, you will see something like the following printed on your screen (some sections may change depending on the PAML version you have installed on your cluster!):

```text
running baseml tmp0001.ctl
BASEML in paml version 4.9h, March 2018
ns = 218        ls = 9869
Reading sequences, sequential format..
Reading seq #218: Tri_lon_Branchiopoda
Sequences read..

9869 site patterns read, 22384 sites
Counting frequencies..

   283836 bytes for distance
 68530336 bytes for conP
   315808 bytes for fhK
  8000000 bytes for space
```

As soon as you see the last line, you will see that various `tmp000X*` files will have been created, and hence you can stop this run by typing `ctrl+C` on the terminal that you have used to run such command. Once you have done this, you can check that the control file you will later need has been created:

```sh
# Run from the `spiders_dating/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
grep 'seqfile' */*/tmp0001.ctl | wc -l # You should get as many datasets as you have, in this case 3
```

Note that, when we ran the commands above, we were not interested in running `BASEML` or `MCMCtree`. We just wanted to execute `MCMCtree` with option `usedata = 3` so that it generates the `tmp000*` files that `BASEML` will later need to estimate the branch lengths, the gradient, and the Hessian. We do this analysis in two steps given that there are restrictions in the HPC we are using that do not allow us to run `BASEML` + `MCMCtree` in one unique job within a reasonable amount of time. In addition, we want to modify some settings in the control file that is automatically generated when enabling `usedata = 3` so that they match what we want to do for our inference. In a nutshell, this is what you will be doing:

1. Run `MCMCtree` to generate the `tmp000*` files.
2. Modify the `tmp0001.ctl` file according to the settings we want to enable to analyse our dataset with `BASEML`.
3. Run `BASEML` using the `tmp000*` files so that it estimates the branch lengths, the gradient, and the Hessian and saves them in a file called `rst2`.
4. Generate the final `in.BV` file for our dataset, which will be later used by `MCMCtree`.

Once all `tmp000*` files are generated for all alignments, we need to make sure that the correct evolutionary model has been enabled (i.e., `model = 4`, `ncatG=4` for HKY85+G4) and that option `method = 1` is enabled, which will speed up the computation of the Hessian and the gradient. We can run the next code snippet to very that the four requirements aforementioned are met:

```sh
# Run from the `spiders_dating/Hessian` dir on your local
# PC. Please change directories until
# you are there. Then, run the following
# commands.
sed -i 's/method\ \=\ 0/method\ \=\ 1/' */*/tmp0001.ctl
grep 'method = 1' */*/tmp0001.ctl | wc -l # You should get as many as datasets you have
grep 'alpha' */*/tmp0001.ctl   # You should see `fix_alpha = 0` and `alpha = 0.5`
grep 'ncatG' */*/tmp0001.ctl   # You should see `ncatG = 4`
grep 'model' */*/tmp0001.ctl   # You should see `model = 3` (i.e., empirical+F model)
```

### Executing `BASEML`

We can now run `BASEML` given that we have the control file ready as well as all the required input files!

We have created a template bash script with flags (i.e., see script `pipeline_Hessian_BASEML_template.sh` in the [`scripts` directory](01_PAML/00_Hessian/scripts/pipeline_Hessian_BASEML_template.sh)), which will be replaced with the appropriate values by another bash script (`generate_job_BASEML.sh`, also saved in the [`scripts` directory](01_PAML/00_Hessian/scripts/generate_job_BASEML.sh)). Please note that the second bash script will edit the template bash script according to the data alignment/s that will be analysed. We had already transferred these scripts to the HPC server when setting up our file structure. Therefore, we just need to execute the following code snippet there:

```sh
# Run from `spiders_dating` dir on your HPC. Please change directories until
# you are there. Then, run the following
# commands.
home_dir=$( pwd )
cd scripts
chmod 775 *sh
num_aln=3
# Arg1: Number of alignments
# Arg2: Path to the pipeline directory
# Arg3: Name of the working directory (i.e., `spiders_dating` in this analysis)
# Arg4: Name of the executable file for BASEML. E.g., `baseml4.10.7`, `baseml`, etc.
# Arg5: Boolean, PAML exported to the path? `Y` of `N`
# Arg6: Requested RAM. E.g., `2G`
./generate_job_BASEML.sh $num_aln $home_dir/pipelines_Hessian spiders_dating baseml4.10.7 N "2G"
```

Next, we will go to the `pipelines_Hessian` directory and run the script that will have been generated using the commands above:

```sh
# Run from `spiders_dating/pipelines_Hessian` dir on your HPC.
# Please change directories until
# you are there. Then, run the following
# commands.
#
# If you list the content of this directory,
# you will see the pipeline you will need 
# to execute in a bash script called
# `pipeline_Hessian.sh`
ll *
# Now, execute this bash script
chmod 775 *sh
qsub pipeline_Hessian.sh
```

Once `BASEML` finishes, we are ready to generate the `in.BV` file that we will later use when running `MCMCtree` to approximate the likelihood calculation:

```sh
# Run from dir `spiders_dating/Hessian/` dir on your HPC
# Please change directories until
# you are there. Then, run the following
# commands.
num_aln=3
for i in `seq 1 $num_aln`
do
printf "\nGenerating in.BV files for dir "$i" ... ...\n\n"
cp $i/rst2 $i/in.BV
done
```

Next, we can transfer the output generated by `BASEML` to our HPC so that we can keep a backup:

```sh
# Run from `00_BASEML_conc` in your PC
mkdir out_BASEML
cd out_BASEML
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/spiders_dating/Hessian .
rsync -avz --copy-links <uname>@<server>:<path_to_your_wd_in_HPC>/spiders_dating/pipelines_Hessian .
# Remove unnecessary empty output files
rm pipelines_Hessian/*sh.o*
```

We can now proceed to timetree inference with `MCMCtree` using the concatenated dataset while enabling cross-bracing across all possible mirrored nodes, regardless a fossil calibrations is constraining their age. [You can click this link to move to the next `README.md` file](../01_MCMCtree/README.md)!
