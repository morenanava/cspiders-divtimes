# Data formatting

Before proceeding with timetree inference, we need to make sure that:

1. The alignment file is in PHYLIP format and easy to read (i.e., ideally one sequence per line).
2. The tree file is in Newick format.

## Alignment files

If you open [the alignment file](00_raw_data/alignment/abce_94_markers_concat_filtered.fasta), you will see that each aligned sequence is already in a unique line. In addition, there are 6 spaces between the sequence ID and the sequence itself, which is the format we need for subsequent formatting steps (see below).

Nevertheless, the sequence IDs are very long, which is not ideal and may lead to errors! We will parse these alignments and generate the equivalent FASTA-formatted files to later fix the sequence tags. For the former, we will use our [PERL in-house script `PHYLtoFASTA.pl`](../src/FASTAtoPHYL.pl), which requires a specific separation (i.e., 6 spaces) between the sequence ID and the sequences:

```sh
# Run the next commands from the 
# `00_data_formatting` directory
cd 00_raw_data/alignment
name=`ls *phy`
sed -i 's/\r//g' *phy
for i in $name
do
txt_name=$( echo $i | sed 's/..*_//' | sed 's/\.phy//' )
printf "Converting "$i" into a FASTA-formatted file\n"
../../../src/PHYLtoFASTA.pl $i
mv log_lenseq.txt log_lenseq_$txt_name".txt"
done
```

Now, we have written an [in-house R script](scripts/Fix_seq_ids.R) to read the FASTA-formatted alignments and the tree files in order to replace the long IDs with shorter ones. Once you run the script, you will see both new FASTA sequence files which name ends with "newIDs.fasta" under the [`alignment` directory](00_raw_data/alignment) and Newick-formatted tree files under the [`trees` directory](00_raw_data/trees) following the same file name. In addition, you can find tsv files under the newly created `out_RData` directory with listing the old tags for each taxon, the various elements comprising such tag, and the newly included tag for `PAML` analyses.

Then, we just need to reconvert the newly FASTA-formatted files into PHYLIP format and save our input files!

```sh
# You should still be inside `00_raw_data/alignment`
# If not, please move to this directory and run the
# following commands
#
# First, fix potential problems with out R files
sed -i 's/\r//g' *newIDs.fasta
# Now, generate PHY format
for aln_name in *newIDs.fasta
do
a_noext=$( echo $aln_name | sed 's/\_newIDs\.fasta//' )
num=$( grep '>' $aln_name | wc -l )
len=$( sed -n '2,2p' $aln_name | sed 's/\r//' | sed 's/\n//' | wc --m )
perl ../../../src/FASTAtoPHYL.pl $aln_name $num $len 
# Create a directory for input data for `MCMCtree`
if [ ! -d "../../01_inp_data" ]
then
mkdir ../../01_inp_data
fi
mv $a_noext"_newIDs.phy" ../../01_inp_data/$a_noext".phy"
rm log_lenseq.txt
done
```

You will now see a new directory called `01_inp_data` inside the directory [`00_data_formatting`](README.md). If you navigate to this newly created `01_inp_data` directory, you will find the alignments in PHYLIP format with the new sequence IDs.

The alignments are now in the correct format, so we can start to parse the tree files!

## Tree files

After using our [R in-house script](scripts/Fix_seq_ids.R) when formatting our alignments, we also formatted our tree files: all the files which name ends with `newIDs.tree` saved inside the [`trees` directory](00_raw_data/trees/) are in Newick format, contain branch lengths, and have the new sequence IDs!

The first thing that we need to do is generating tree file with only the tree topology (i.e., no branch lengths) so that we can run `CODEML` to calculate the Hessian, the gradient, and the branch lengths prior to timetree inference -- they are needed to approximate the likelihood calculation!

```sh
# Run from `00_raw_data/trees`, newly created directory
cp tree_noeurycyde_calibnames_newIDs.tree ../../01_inp_data/tree_noeurycyde_uncalib.tree
cp tree_noUCEs_calibnames_newIDs.tree ../../01_inp_data/tree_noUCEs_uncalib.tree
cp tree_supaln_calibnames_newIDs.tree ../../01_inp_data/tree_supaln_uncalib.tree
sed -i 's/e\-//g' ../../01_inp_data/*uncalib.tree
sed -i 's/:[0-9]*\.[0-9]*//g' ../../01_inp_data/*uncalib.tree
sed -i 's/:[0-9]*//g' ../../01_inp_data/*uncalib.tree
# Add headers
sed -i '1s/^/218 1\n/' ../../01_inp_data/*uncalib.tree
sed -i '1s/^/218 1\n/' LUCAdup_allcb_calibnames.tree
sed -i '1s/^/218 1\n/' LUCAdup_allcb_calibnames.tree
```

Now, we can calibrate these topologies following the calibration files we have designed. In a nutshell, the tasks carried out by our [R in-house script](scripts/Include_calibrations.R) are the following:

* Read the input file with calibration information (i.e., [`Calibs_*txt` files inside `calibs`, one per alignment](00_raw_data/calibs/)), and the corresponding input tree files with the tree topology that shall be fixed when running `PAML` programs (i.e., the uncalibrated tree files we have just generated).
* Find the corresponding nodes for the MRCAs for each calibrated node so that the tag that will be used to identify the calibration for such node (i.e., first column in the calibration files) is incorporated as a node label.
* Identify nodes that are assigned more than one calibration (i.e., potential duplicates in the calibration file) and decide what to do with them. Then, repeat the previous step with the filtered calibration file, if generated, and output the resulting tree topology with node labels.
* Read the previously output file as a character vector (not as a `phylo` object!) to replace the node labels with the corresponding node calibrations in `MCMCtree` notation. Columns 4-7 in the calibration files are used for this purpose.

Once you have run our [R in-house script `Include_calibrations.R`](scripts/Include_calibrations.R), you will see that the calibrated tree files with nodes labelled following `MCMCtree` notation are output in the [`01_inp_data`](01_inp_data) directory (i.e., check all files which name ends with `calib_MCMCtree.tree`). In addition, files ending with `fordisplay_calib_MCMCtree.tree` will have been output in directory [`00_raw_data/trees`](00_raw_data/trees/), which can be opened with graphical viewers such as `FigTree` or `TreeViewer` ([Bianchini and Sánchez-Baracaldo, 2024](https://onlinelibrary.wiley.com/doi/10.1002/ece3.10873)) to see which nodes have been calibrated and how.

Now, we can move onto the next step: [we can calculate the Hessian, the branch lengths, and the gradient with `CODEML`!](../01_PAML/00_CODEML/README.md)
