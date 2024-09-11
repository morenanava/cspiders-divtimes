# A timescale to sea spiders' evolution -- timetree inference analyses

In this repository, you will find a very detailed tutorial with all the steps you need to follow to reproduce the results for the timetree inference analyses we carried out as part of our study: **A timescale to sea spiders' evolution**.

To get started, you can clone this repository on your PC and follow all the guidelines given in the various `README.md` files that you shall find inside each directory so that you can go through all the steps we carried out for timetree inference. A summary of this workflow is given below:

* Parsing and formatting the input data required to run `PAML` programs `BASEML` and `MCMCtree` ([Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17483113/)): sequence files, calibrated tree files (more details below), and control files.
* Inferring the mean evolutionary rate to specify a sensible rate prior.
* Running `PAML` programs for timetree inference:
  * Using various in-house pipelines to set up the working environment, the file structure, and the control files required to run `PAML` programs.
  * Running `CODEML` to calculate the branch lengths, the gradient, and the Hessian; required by `MCMCtree` to enable the approximate likelihood calculation [dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613).
  * Running `MCMCtree` with the approximate likelihood calculation enabled for timetree inference. We assessed the impact that **including UCEs in the molecular alignment** and **constraining a node age using a fossil specimen still under debate** could have on timetree inference under both the Geometric Brownian motion (**GBM**, [Thorne et al., 1998](https://pubmed.ncbi.nlm.nih.gov/9866200/), [Yang and Rannala, 2006](https://pubmed.ncbi.nlm.nih.gov/16177230/)) and the independent-rates log-normal (**ILN**, [Rannala and Yang, 2007](https://pubmed.ncbi.nlm.nih.gov/17558967/), [Lemey et al., 2010](https://pubmed.ncbi.nlm.nih.gov/20203288/)) **relaxed-clock models**. The following tags are used to name files, directories, or specific analyses throughout this tutorial depending on the sequence alignment used and the calibration strategy chosen:
  * **`noeurycyde`**:
    * The sequence alignment includes the mitochondrial, 18S, and UCE sequences (i.e., what we label as "Matrix 1" in our manuscript).
    * This calibration strategy does not use *?Eurycyde golem* to constrain the age of node Ascorhynchidae. A total of four node ages are constrained: Arthropoda (root), Chelicerata, Colossendeidae, Phoxichilidioidea.
  * **`noUCEs`**:
    * The sequence alignment includes only the mitochondrial and 18S sequences (i.e., what we label as "Matrix 2" in our manuscript, UCEs were excluded).
    * The calibration strategy uses *?Eurycyde golem* to constrain the age of node Ascorhynchidae. A total of five node ages are constrained: Arthropoda (root), Chelicerata, Colossendeidae, Phoxichilidioidea, and Ascorhynchidae.
  * **`supaln`**:
    * The sequence alignment includes the mitochondrial, 18S, and UCE sequences (i.e., what we label as "Matrix 1" in our manuscript).
    * The calibration strategy uses *?Eurycyde golem* to constrain the age of node Ascorhynchidae. A total of five node ages are constrained: Arthropoda (root), Chelicerata, Colossendeidae, Phoxichilidioidea, and Ascorhynchidae.
* Running MCMC diagnostics for all the chains under each analysis.

To make it easier for you to navigate this GitHub repository, below you can find a summary of the content you shall find inside each directory.

## [Data formatting](00_data_formatting/README.md)

* **Molecular alignment**: a total of 13 mitochondrial protein-coding genes, 18S rRNA (18S) nuclear genes, and 98 Ultra Conserved Elements (UCEs) were retrieved from the sequences made available by [Ballesteros et al. (2021)](https://academic.oup.com/mbe/article/38/2/686/5904272) and [Sabroux et al. (2023a)](https://www.sciencedirect.com/science/article/abs/pii/S105579032300026X). The final dataset entails 198 pycnogonid taxa across the 11 families currently accepted in Pantopoda ([Bamber et al. 2024](https://mapress.com/zt/article/view/zootaxa.1668.1.15)). Outgroups include thirteen non-pycnogonid chelicerates, four pancrustaceans, and three myriapods. Scorpions, non-Mesothelae spiders, and Trombidoformes mites were excluded from our analyses in order to prevent the potential emergence of tree reconstruction artefacts caused by convergent inversions of base composition bias in their mitochondrial genomes ([Hassanin et al. 2005](https://pubmed.ncbi.nlm.nih.gov/16021696/), [Arabi et al. 2012](https://pubmed.ncbi.nlm.nih.gov/22362465/), [Sabroux et al. 2023a](https://www.sciencedirect.com/science/article/abs/pii/S105579032300026X)). Following [Rota-Stabelli and Telford (2008)](https://pubmed.ncbi.nlm.nih.gov/18501642/), we did not include Pseudoscorpiones nor Mesostigmata because these taxa yield long branches in mitogenomic trees. Given that UCEs were only available for a subset of the species present in our dataset, we decided to account for the potential impact that a high rate of missing data in UCEs (95.70% on average per specimen) could have on phylogeny and timetree inference analyses. We therefore built two molecular alignments:
  * **"Matrix 1**: mitochondrial, 18S, and UCE sequences were included. `MAFFT-linsi` v.7.471 ([Katoh & Standley 2013](https://pubmed.ncbi.nlm.nih.gov/23329690/)) and trimmed with `trimAl` v.1.2 (`Gappyout` option; [Capella-Gutiérrez et al. 2009](https://pubmed.ncbi.nlm.nih.gov/19505945/)). The size of “Matrix 1” was reduced by approximately 38.4% (from 36,330 base pairs, to 22,384 bp): mitochondrial (5,579 bp), 18S (926 bp), and UCE (15,879 bp). The tags used throughout this tutorial that indicate that this alignment was used are `supaln` and `noeurycyde`.
  * **"Matrix 2**: mitochondrial and 18S sequences were included. `MAFFT-linsi` v.7.471 ([Katoh & Standley 2013](https://pubmed.ncbi.nlm.nih.gov/23329690/)) and trimmed with `trimAl` v.1.2 (`Gappyout` option; [Capella-Gutiérrez et al. 2009](https://pubmed.ncbi.nlm.nih.gov/19505945/)). The size of “Matrix 2” was reduced by approximately 62.4% (from 17,284 bp, to 6,505 bp): mitochondrial (5,579 bp) and 18S (926 bp). The tag used throughout this tutorial that indicate that this alignment was used is `noUCEs`.
Below, you can find shortcuts to the unformatted alignment files and the final input sequence files:
  * **Unformatted alignments**: directory [`00_raw_data/alignment`](00_data_formatting/00_raw_data/alignment/), inside [`00_data_formatting`](00_data_formatting), contains the alignments in FASTA and PHYLIP format that had a very long tag ID that would cause downstream issues if used with `PAML` programs. Those file names ending with `newIDs` were the ones later converted to PHYLIP format as input sequence files (see [`README.md`](00_data_formatting/README.md#alignment-files) and our in-house R script [`Fix_seq_ids.R`](00_data_formatting/scripts/Fix_seq_ids.R) to fix the tag IDs).
  * **Formatted alignments**: directory [`01_inp_data`](00_data_formatting/01_inp_data/) contains all input sequence files that were used in all timetree inference analyses. Despite "Matrix 1" was used in two analyses (i.e., when either including or removing the age constraint on node Ascorhynchidae), we decided to generate two different sequence files to distinguish the analysis in which such a file was being used.
* **Phylogeny**: we used `IQ-TREE` v2.1.3 ([Minh et al. 2020](https://pubmed.ncbi.nlm.nih.gov/32011700/)) using the ultrafast bootstrap approximation method ([Hoang et al. 2018](https://academic.oup.com/mbe/article/35/2/518/4565479); option `-B 1000` to generate 1,000 replicates) to infer the tree topology that we fixed in all timetree inference analyses. Each superalignment (i.e., "Matrix 1" and "Matrix 2") was partitioned by genes during phylogeny inference to account for the varying rates of evolution of each gene of our dataset (option `-p` enabled). `ModelFinder` (option `-m MFP`) was used to find the best-fitting substitution models ([Kalyaanamoorthy et al. 2017](https://www.nature.com/articles/nmeth.4285)) under which phylogeny inference took place. Given that we wanted to test two different calibration strategies with "Matrix 1", we decided to have two tree input files and two alignment input files despite the latter were the same (i.e., only the file name changed to make it easier for our in-house scripts to semi-automate the analyses; see [`tree_noeurycyde_calibnames_newIDs.tree`](00_data_formatting/00_raw_data/trees/tree_noeurycyde_calibnames_newIDs.tree) and [`tree_supaln_calibnames_newIDs.tree`](00_data_formatting/00_raw_data/trees/tree_supaln_calibnames_newIDs.tree)). We then used our R in-house script [`Include_calibrations.R`](00_data_formatting/scripts/Include_calibrations.R), together with our input files with information about the node age constraints (see directory [`calibs`](00_data_formatting/00_raw_data/calibs/)), to obtain three calibrated tree files:
  * **[`noeurycyde_calib_MCMCtree.tree`](00_data_formatting/01_inp_data/noeurycyde_calib_MCMCtree.tree)**: calibrated tree file with four node age constraints (i.e., *?Eurycyde golem* was not used). Timetree inference analyses using this input tree file require the usage of "Matrix 1" in file [`aln_noeurycyde.phy`](00_data_formatting/01_inp_data/aln_noeurycyde.phy).
  * **[`noUCEs_calib_MCMCtree.tree`](00_data_formatting/01_inp_data/noUCEs_calib_MCMCtree.tree)**: calibrated tree file with five node age constraints (i.e., *?Eurycyde golem* was used). Timetree inference analyses using this input tree file require the usage of "Matrix 2" in file [`aln_noUCEs.phy`](00_data_formatting/01_inp_data/aln_noUCEs.phy).
  * **[`supaln_calib_MCMCtree.tree`](00_data_formatting/01_inp_data/supaln_calib_MCMCtree.tree)**: calibrated tree file with five node age constraints (i.e., *?Eurycyde golem* was used). Timetree inference analyses using this input tree file require the usage of "Matrix 1" in file [`aln_supaln.phy`](00_data_formatting/01_inp_data/aln_supaln.phy).
We then generated the uncalibrated tree file to be used with `BASEML` to infer the branch lengths, the gradient, and the Hessian required to approximate the likelihood calculation. The tree file names that end with `_uncalib.tree` in [`01_inp_data`](00_data_formatting/01_inp_data/) follow the same format as described above for the uncalibrated tree files.
* **Calibrations**: we created one file for each type of analysis as aforementioned: [`Calibs_noeurycyde.txt`](00_data_formatting/00_raw_data/calibs/Calibs_noeurycyde.txt) and [`Calibs_noUCEs.txt`](00_data_formatting/00_raw_data/calibs/Calibs_noUCEs.txt), [`Calibs_supaln.txt`](00_data_formatting/00_raw_data/calibs/Calibs_supaln.txt). These files were formatted so that our [R in-house script `Include_calibrations_MCMCtree.R`](00_data_formatting/scripts/Include_calibrations_MCMCtree.R) could constrain the age of the relevant nodes. There are intermediate files generated by this script inside directory [`trees`](00_data_formatting/00_raw_data/trees), but only the calibrated tree files saved in directory `01_inp_data` are relevant (see links above).

> [!IMPORTANT]
> Please read the [`README.md` file](00_data_formatting/README.md) inside the [`00_data_formatting` directory](00_data_formatting), which has all the information you need to understand how we formatted our data, including code snippets that you can run while following the step-by-step procedure and links to our R in-house scripts.

## [`PAML` analyses](01_PAML)

Below, you can find a short description of what you can find in each directory:

* [**`00_CODEML`**](01_PAML/00_CODEML/README.md): guidelines to infer the branch lengths, the gradient, and the Hessian with `CODEML`.
* [**`01_MCMCtree`**](01_PAML/01_MCMCtree/README.md): guidelines to run `MCMCtree` by enabling the approximate likelihood calculation ([dos Reis and Yang, 2011](https://academic.oup.com/mbe/article/28/7/2161/1051613)) under the **GBM and ILN relaxed-clock models** and when assessing the **impact of different calibration strategies and the inclusion of UCEs**. The directories with name `GBM` correspond to analyses under the autocorrelated-rates model and those with `ILN` under the independent-rates log-normal rates model. Inside these two directories, you will have directories `1`, `2`, and `3`:
  * Directory `1` contains all the relevant files for those analyses when the age of node Ascorhynchidae was **NOT** constrained and "Matrix 1" (UCEs **WERE** included) was used (i.e., file names will include `noeruycyde`).
  * Directory `2` will have those files generated and/or used when the age of node Ascorhynchidae **WAS** constrained and "Matrix 2" (UCEs were **NOT** included) (i.e., file names will include `noUCEs`).
  * Directory `3` will include those files generated and/or used when the age of node Ascorhynchidae **WAS** constrained and "Matrix 1" (UCEs **WERE** included) (i.e., file names will include `supaln`).
A total of 16 chains were run when sampling from the posterior, and so you shall find directories from `1` to `16` inside directories [`GBM/[1-3]`](01_PAML/01_MCMCtree/sum_analyses/01_posterior/GBM) and [`ILN/[1-3]`](01_PAML/01_MCMCtree/sum_analyses/01_posterior/ILN). When sampling from the prior, a total of 6 chains were run, and so you shall find directories from `1` to `6` inside directory [`CLK/[1-3]`](01_PAML/01_MCMCtree/sum_analyses/00_prior/CLK). All the plots generated as part of the MCMC diagnostics can be found under directory [`plots`](01_PAML/01_MCMCtree/plots/); tables with convergence, efficiency, and sampling diagnostics under directory [`out_RData`](01_PAML/01_MCMCtree/out_RData/); and "summary directories" [`sum_files_post`](01_PAML/01_MCMCtree/sum_files_post/) and [`sum_files_prior`](01_PAML/01_MCMCtree/sum_files_prior/) have the most relevant files generated during MCMC diagnostics when sampling from the posterior and the prior, respectively. The output tab-separated files that were used in subsequent steps of this study are the following:
  * [**Summary of estimated divergence times | `noeurycyde + GBM`**](01_PAML/01_MCMCtree/sum_files_post/noeurycyde_GBM_all_mean_est.tsv): this file has the estimated mean posterior divergence times (second column), 2.5% quantile (third column), and 97.5% quantile (fourth column) for those analyses ran under the **GBM relaxed-clock model**, when the age of node "Ascorhynchidae" was **not** constrained, and when **"Matrix 1"** (UCEs included) was used.
  * [**Summary of estimated divergence times | `noeurycyde + ILN`**](01_PAML/01_MCMCtree/sum_files_post/noeurycyde_ILN_all_mean_est.tsv): this file has the estimated mean posterior divergence times (second column), 2.5% quantile (third column), and 97.5% quantile (fourth column) for those analyses ran under the **ILN relaxed-clock model**, when the age of node "Ascorhynchidae" was not constrained, and when **"Matrix 1"** (UCEs included) was used.
  * [**Summary of estimated divergence times | `noUCEs + GBM`**](01_PAML/01_MCMCtree/sum_files_post/noUCEs_GBM_all_mean_est.tsv): this file has the estimated mean posterior divergence times (second column), 2.5% quantile (third column), and 97.5% quantile (fourth column) for those analyses ran under the **GBM relaxed-clock model**, when the age of node "Ascorhynchidae" **was** constrained, and when **"Matrix 2"** (UCEs excluded) was used.
  * [**Summary of estimated divergence times | `noUCEs + ILN`**](01_PAML/01_MCMCtree/sum_files_post/noUCEs_ILN_all_mean_est.tsv): this file has the estimated mean posterior divergence times (second column), 2.5% quantile (third column), and 97.5% quantile (fourth column) for those analyses ran under the **ILN relaxed-clock model**, when the age of node "Ascorhynchidae" **was** constrained, and when **"Matrix 2"** (UCEs excluded) was used.
  * [**Summary of estimated divergence times | `supaln + GBM`**](01_PAML/01_MCMCtree/sum_files_post/supaln_GBM_all_mean_est.tsv): this file has the estimated mean posterior divergence times (second column), 2.5% quantile (third column), and 97.5% quantile (fourth column) for those analyses ran under the **GBM relaxed-clock model**, when the age of node "Ascorhynchidae" **was** constrained, and when **"Matrix 1"** (UCEs included) was used.
  * [**Summary of estimated divergence times | `supaln + ILN`**](01_PAML/01_MCMCtree/sum_files_post/supaln_ILN_all_mean_est.tsv): this file has the estimated mean posterior divergence times (second column), 2.5% quantile (third column), and 97.5% quantile (fourth column) for those analyses ran under the **ILN relaxed-clock model**, when the age of node "Ascorhynchidae" **was** constrained, , and when **"Matrix 1"** (UCEs uncluded) was used.

Please note that each directory has all the data and scripts that you will need to reproduce our analyses if you clone this repository on your PC. Nevertheless, note that the `MCMCtree` output files with all the collected samples during the MCMC that we used to summarise our results are too large to be stored here. If you want to use our results (e.g., compare to the results you get when you re-run our analyses, reproduce our figures or plots, etc.), please download [the GitHub repository from our data archive](.). Depending on how you want to use our step-by-step guidelines, you may do one of the following:

* If you want to run the tutorials under the `01_PAML/01_MCMCtree` to reproduce the results we report in our manuscript, then you will not initially need our output files: you will run `PAML` programs following our guidelines to generate your own results, which you can analyse with our scripts. Please bear in mind that the your time estimates may vary (i.e., decimals) because you will be using different seed numbers: you will be carrying out independent runs with `MCMCtree`, and such variation is expected (unless the same seed number is used!).
* If you want to either (i) reproduce our figures or plots using our own results or (ii) you have re-run our analysis and want to compare your results to ours, please make sure that you do the following:
  * Download the compressed file with all our timetree inference results from [our data archive)](.). Please note that we may have changed the content of some of the `README.md` files and/or other files since we generated our archive, and so we ask that you use such archive to access the files that you may not find in this repository (due to size limitation) if you need them for the reasons aforementioned. **Please keep an eye on this GitHub repository and our releases to track all the changes/updates since the publication of our paper**.
  * Once downloaded, please decompress the `cspiders-divtimes.zip` file and find the `00_prior` and `01_posterior` directories inside each `01_PAML/01_MCMCtree/sum_analyses` directory. You will see that the `00_prior` and `01_posterior` directories are already inside the `01_PAML/01_MCMCtree` directories that we provide you with in this GitHub repository. Nevertheless, they are all missing the `mcmc.txt` files that you shall find inside the downloaded `00_prior/CLK/*/` and `01_posterior/[GBM|ILN]/*/` directories, which are the output files generated by `MCMCtree` with all the samples collected for each parameter of interest. Our in-house scripts use such files to summarise the results. We recommend that you do the following:
    * If you just want to recreate figures/tables/plots we have included in our paper with our results, please navigate to `01_PAML/01_MCMCtree` in the decompressed file you should have already downloaded, locate the `00_prior` and `01_posterior` directories, and save them in the corresponding `sum_analyses` directories in the file structure you are following as per this repository. Make sure that you save them in the same `01_PAML/01_MCMCtre` directory following the file structure! Then, you can run our scripts to generate the files you wanted.
    * If you want to compare your results to ours, please change the name of your `sum_analyses` directory with your own results so that you do not overwrite your files. Then, you can locate our `sum_analyses` directories in the downloaded archive and place them in the corresponding location in your file structure. You can then re-run our scripts with both your results and ours if you want to also test that you can recreate the same results that we report in our paper.

In a nutshell, if you follow the instructions and run the code snippets and in-house scripts as detailed in each `README.md` file under directories `00_CODEML` and `01_MCMCtree`, you will be able to reproduce all our results! If you have re-run these analyses, then you shall be able to compare your results to ours!

> [!IMPORTANT]
> Note that all the directories that you shall find inside the [`01_PAML` directory](01_PAML) have their own `README.md` file. To this end, you can navigate to each of these directories to read more about how we ran `CODEML` and `MCMCtree` under different settings (e.g., data type, calibration strategy, relaxed-clock model).

## [Range Through Time plot](02_RTT_plots/README.md)

The in-house R script [`Create_RTT_plot.R`](scripts/Create_RTT_plot.R) is used to generate the range through time plot that shows the distribution of the sea spider fossil record. The function used to create the plot modifies an original function within the R package `paleoverse` so that the range plot can be ordered by groups above species. The in-house function can be found in [`Functions_RTTplot.R](scripts/Functions_RTTplot.R). More details on the main [`README.md` file](02_RTT_plots/README.md).

## [In-house scripts](src)

The main scripts that you call when you follow this tutorial are detailed below:

* [Functions.R](src/Functions.R): this script has all the main in-house function we have written to analyse timetree inference results generated with `MCMCtree`. When running the MCMC diagnostics, you will be using those!
* [FASTAtoPHYL.pl](src/FASTAtoPHYL.pl): this script formats FASTA sequence files into PHYLIP sequence files.
* [PHYLtoFASTA.pl](src/PHYLtoFASTA.pl): this script formats PHYLIP sequence files into FASTA sequence files.
* [one_line_fasta.pl](src/one_line_fasta.pl): this script makes sure that there is one line for each sequence in FASTA sequence files.

If you open the scripts aforementioned, you will find all the details regarding the arguments they need and how to run them (which is also referred to when they are called when summarising our results and/or generating output files such as plots or tables).

## What do you need to run our analyses?

### Software

Before you go through this step-by-step tutorial to reproduce our results, please make sure you have the following software installed on your PCs/HPCs:

* **`PAML`**: we used `PAML` v4.10.7, available from the [`PAML` GitHub repository](https://github.com/abacus-gene/paml). Please download [the latest release on the PAML GitHub repository which, as of January 2024, corresponds to the `PAML` version we used: v4.10.7](https://github.com/abacus-gene/paml/releases/tag/4.10.7). You can use either the pre-compiled binaries for your OS or compile the source code. Note that we used the Linux version to run these analyses on an HPC. Numerical differences may occur depending on the OS where you run this software.

> [!NOTE]
>
> **Linux users**: you may need to install the `intel` compilers before you run `make -f Makefile`. Please visit this link to download the [Intel oneAPI Base Toolkit](https://www.intel.com/content/www/us/en/developer/articles/tool/oneapi-standalone-components.html#dpcpp-cpp). Thank you, [@ERRMoody](https://github.com/ERRMoody), for pointing this out!

* **`R`** and **`RStudio`**: please download [R](https://cran.r-project.org/) (i.e., select either `Download R for Linux`, `Download R for macOS`, or `Download R for Windows` depending on your OS; then follow the installation instructions) and [RStudio](https://posit.co/download/rstudio-desktop/) as you will need them to run various of our in-house R scripts. The packages you will need to install work with R versions that are either newer than or equal to v4.1.2 (we used v4.1.2). If you are a Windows user, please make sure that you have the correct version of `RTools` installed, which will allow you to install packages from the source code if required. For instance, if you have R v4.1.2, then installing `RTools4.0` shall be fine. If you have another R version installed on your PC, please check whether you need to install `RTools 4.2` or `RTools 4.3`. For more information on which version you should download, [please go to the CRAN website by following this link and download the version you need](https://cran.r-project.org/bin/windows/Rtools/).

    Before you proceed, however, please make sure that you install the following packages too:

    ```R
    # Run from the R console in RStudio
    # Check that you have at least R v4.1.2
    version$version.string
    # Now, install the packages you will need to run
    # our in-house R scripts
    # Note that it may take a while to install them all
    install.packages( c('rstudioapi', 'ape', 'phytools', 'sn', 'stringr', 'rstan'), dep = TRUE )
    ## NOTE: If you are a Windows user and see the message "Do you want to install from sources the 
    ## packages which need compilarion?", please make sure that you have installed the `RTools`
    ## aforementioned. Otherwise, you will get errors during the installation.
    ## The versions we used for each of the packages aforementioned are the following:
    ##   rstudioapi: v0.14
    ##   ape: v5.7.1
    ##   phytools: v1.5.1
    ##   sn: v2.1.1
    ##   stringr: v1.5.0
    ##   rstan: v2.21.7
    ```
  
> [!NOTE]
>
> **Linux users**: you may experience some problems when you try to execute the `install.packages()` command that you see in the code snippet above. If that is the case, please install each package separately (e.g., `install.packages( 'rstudioapi' )` to install `rstudioapi`, etc.). If you experience problems with `stringr`, please follow the [suggestions given in this StackOverflow page](https://stackoverflow.com/questions/38987157/libicu-and-stringi-on-fedora-24-causing-r-headaches/39411793#39411793) -- despite the question being addressed for Fedora, the solution suggested also works for Ubuntu. Thank you, [@ERRMoody](https://github.com/ERRMoody), for pointing this out!

* **`TreeViewer`**: you can use `TreeViewer` ([Bianchini and Sánchez-Baracaldo, 2024](https://onlinelibrary.wiley.com/doi/10.1002/ece3.10873)) as a graphical interface with which you can display and highly customise the format of the timetrees we have generated and include in this repository. You may want to read [the `TreeViewer` documentation](https://github.com/arklumpus/TreeViewer/wiki) to learn more about which modules you can use and how you can improve the design of your timetrees (e.g., include/exclude densities, include pictures, play with various colours and shapes, etc.). You can [download `TreeViewer` by following this link](https://treeviewer.org/).

* **`FigTree`**: alternatively, you can use `FigTree` to display the timetrees we have generated. While not as customisable as `TreeViewer`, you can decide what you want to be displayed on the graphical interface by selecting the buttons and options that enable a specific design. You can [download the latest pre-compiled binaries, `FigTree v1.4.4`, from the `FigTree` GitHub repository](https://github.com/rambaut/figtree/releases).

* **`Tracer`**: you can use this graphical interface to visually assess the MCMCs you have run during your analyses (e.g., chain efficiency, chain convergence, autocorrelation, etc.). You can [download the latest pre-compiled binaries, `Tracer v1.7.2`, from the `Tracer` GitHub repository](https://github.com/beast-dev/tracer/releases/tag/v1.7.2).

### Additional note for Mac users when using the `sed` command

If you are running a UNIX-based system (e.g., Mac users), you will experience some problems with the Linux-based `sed` command that you shall see (i) in the various code snippets that you need to run as part of the tutorial described in this repository and (ii) used in most of our in-house bash scripts. By default, this command is different from Linux-based systems, and hence you will have problems to execute it as intended unless you follow one of these two approaches:

* (EASY): Instead of running the `sed` commands using the format `sed -i 's/PATTERN/REPLACEMENT/'` (i.e., what you shall see in the code snippets), you should include `''` between `-i` and `'s/PATTERN/REPLACEMENT/'`: `sed -i '' 's/PATTERN/REPLACEMENT/'`. Please remember to modify the commands in this tutorial accordingly before you copy and paste them on your terminal!
* (MORE DIFFICULT): You should install `GNU sed` and establish it as the "standard" `sed` command instead of the one you will have by default. At the time of writing, [this post](https://medium.com/@bramblexu/install-gnu-sed-on-mac-os-and-set-it-as-default-7c17ef1b8f64) is available and has a detailed explanation that you could follow to install `GNU sed`. Nevertheless, there are many other tutorials out there that you can follow to achieve the same goal. If you decide to follow this approach, you will not need to change the commands in this tutorial!

Once the `sed` incompatibility is sorted out, there should not be any problems with regards to following this tutorial and/or running our in-house bash scripts!

<br>

----

<br>

If you have gone through the previous sections and have a clear understanding of the data we used, the workflow we followed (which you shall follow to reproduce our analyses), and have installed the required software aforementioned... Then you are ready to go!

You can start by taking a look at how we formatted the raw dataset with which we started the timetree inference analyses [by following this link](00_data_formatting/README.md).

Happy timetree inference! :)

## Contact

* **Morena Nava** ([`@morenanava`](https://github.com/morenanava/)).
* **Sandra Álvarez-Carretero** ([`@sabifo4`](https://github.com/sabifo4/)).

If you have any queries with regards to the analyses detailed in this repository, please do not hesitate to <a href="mailto:sandra.ac93@gmail.com">get in touch via email</a>!