#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
# Load package to help find the path to this source file 
library(rstudioapi) 
# Get the path of current open file
path_to_file <- getActiveDocumentContext()$path 
# Get working directory path
wd      <- paste( dirname( path_to_file ), "/", sep = "" )
wd.name <- dirname( path_to_file )
# Set wd. Note that the wd will be the same directory where this 
# script is saved.
setwd( wd )
# Load in-house R script with main function to add node age constraints
source( "../../src/Functions.R" )

#----------------------------#
# DEFINE GLOBAL VARS BY USER #
#----------------------------#
# Path to your input tree with calibrations. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file. You need to include the
# flags within square brackets (e.g., [Mammalia]) and write them on the node
# that is to be calibrated.
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.


# Path to your input text file that allows you to match the flags you have 
# used to constrain node ages with the calibration you 
# want to use in `MCMCtree` notation. The format you need to follow is given
# below:
#
#   - Header.
#   - One row per calibration.
#   - No spaces at all, semi-colon separated.
#   - There are 4 columns:
#       - Name you want to give to the calibrated node (no spaces!).
#       - Name of one of the tips (e.g., tip 1) that leads to MRCA (no spaces!).
#       - Name of the other tip (e.g., tip 2) that leads to MRCA (no spaces!).
#       - Calibration in `MCMCtree` notation (no spaces!). More details on the 
#         `MCMCtree` notation you need to use in the fourth column in the PAML
#         documentation:
#         https://github.com/abacus-gene/paml/blob/master/doc/pamlDOC.pdf
# 
# E.g.: Header and one row in a calibration text file:
#
# ```
# name;tip1;tip2;MCMCtree
# root;sp1;sp2;'B(0.256,1.34,0.025,1e-300)'
# ```
#
# NOTE: Always check that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files. This file needs to be
# already in PHYLIP format. Please follow the same format as used in the 
# example tree file provided.
path_textconv <- c( "../00_raw_data/calibs/Calibs_noeurycyde.txt",
                    "../00_raw_data/calibs/Calibs_noUCEs.txt",
                    "../00_raw_data/calibs/Calibs_supaln.txt"
                    )
all_calibs    <- vector( mode = "list", length(path_textconv) )
names( all_calibs ) <- c( "noeurycyde", "noUCEs", "supaln" )
for( i in 1:length( all_calibs ) ){
  all_calibs[[ i ]] <- read.table( file = path_textconv[i],
                                   stringsAsFactors = FALSE, sep = ";",
                                   blank.lines.skip = TRUE, header = TRUE,
                                   colClasses = rep( "character", 4 ) )
}
# Path to trees
path_trees <- c( "../01_inp_data/tree_noeurycyde_uncalib.tree",
                 "../01_inp_data/tree_noUCEs_uncalib.tree",
                 "../01_inp_data/tree_supaln_uncalib.tree" )
# Run `add_node_const` function to calibrate the tree
# 
# Arguments:
# tt            Phylo, tree object.
# calibrations  Matrix, calibration info, extract from object `all_calibs`. 
#               E.g., `all_calibs[[1]]`
# out_name      Character, name of the dataset analysed.
# out_dir_raw   Character, abs/rel path to the output directory where 
#               calibrated trees that can be visualised by graphical interfaces
#               will be output alongside tmp trees that will be used by
#               this function to generate input data for PAML programs.
# out_dir_inp   Character, abs/rel path to the output directory where inp data
#               that will be used by PAML programs will be saved.
for( i in 1:length( path_trees ) ){
  tt_ape <- ape::read.tree( file = path_trees[i] )
  test   <- add_node_const( tt = tt_ape, calibrations = all_calibs[[ i ]],
                            out_name = names( all_calibs )[i],
                            out_dir_raw = "../00_raw_data/trees/",
                            out_dir_inp = "../01_inp_data/")
  
}

