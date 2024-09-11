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

#-------#
# TASKS #
#-------#
# Get sequence files
fasta_fpath  <- list.files( path = "../00_raw_data/alignment/", 
                            pattern = "fasta", full.names = TRUE )
##> CHECK
del_fpath    <- grep( pattern = "newIDs", x = fasta_fpath )
if( length( del_fpath) > 0 ){
  fasta_fpath <- fasta_fpath[-c(del_fpath)]
}
##> END CHECK
fasta_fnames <- gsub( pattern = "..*\\/|.fasta", replacement = "" ,
                      x = fasta_fpath )
fasta_ff          <- vector( mode = "list", length = length( fasta_fnames ) )
names( fasta_ff ) <- fasta_fnames
# Get tree files
tree_path <- list.files( path = "../00_raw_data/trees/",
                         pattern = "calibnames.tree",
                         full.names = TRUE )
##> CHECK
del_tpath    <- grep( pattern = "newIDs", x = tree_path )
if( length( del_tpath) > 0 ){
  tree_path <- fasta_fpath[-c(del_tpath)]
}
##> END CHECK
# Save files
for( i in 1:length( fasta_ff ) ){
  tmp_ff <- ape::read.dna( file = fasta_fpath[i], format = "fasta" )
  # Get seq IDs, and remove `_concat`
  rows_nraw <- rownames( tmp_ff )
  rows_n <- gsub( x = rows_nraw, pattern = "_concat", replacement = "" ) 
  # Create matrix that will be output with a summary later
  out_mat <- matrix( 0, nrow = length(rows_n), ncol = 6 )
  colnames( out_mat ) <- c( "Genus", "2nd_element", "3rd_element",
                            "4th_element", "PAML_tag", "Old_tag" )
  out_mat[,6] <- rows_nraw
  # Start populating matrix and generating new tag
  p_sep4 <- stringr::str_split_fixed( string = rows_n, pattern = "_", n = 4 )
  p1 <- stringr::str_sub( string = p_sep4[,1], 1, 3 )
  out_mat[,1] <- p_sep4[,1]
  p2 <- stringr::str_sub( string = p_sep4[,2], 1, 3 )
  out_mat[,2] <- p_sep4[,2]
  p3 <- gsub( pattern = "_", replacement = "", x = p_sep4[,3] )
  out_mat[,3] <- p_sep4[,3]
  p4 <- gsub( pattern = "_", replacement = "", x = p_sep4[,4] )
  out_mat[,4] <- p_sep4[,4]
  tmp_tag <- vector( mode = "character", length = dim(p_sep4)[1] )
  count <- 0
  track_visited <- rep( 0, length = dim(p_sep4)[1] )
  for( j in p3 ){
    count <- count + 1
    if( j == "" & p4[count] == "" ){
      tmp_tag[count] <- gsub( x = paste( p1[count], "_", p_sep4[count,2],
                                         sep = "" ),
                              pattern = "\\.", replacement = "" )
      track_visited[count] <- count
      out_mat[count,5] <- tmp_tag[count]
    }
  }
  rm_tv <- which( track_visited == 0 )
  track_visited <- track_visited[-c(rm_tv)]
  left_pos <- c(1:dim(p_sep4)[1])[-track_visited]
  for( j in left_pos ){
    if( p4[j] == "" ){
      tmp_tag[j] <- gsub( x = paste( p1[j], "_", p2[j], "_",
                                         p3[j], sep = "" ),
                              pattern = "\\.", replacement = "" )
      out_mat[j,5] <- tmp_tag[j]
    }else{
      tmp_p3     <- stringr::str_sub( string = p_sep4[j,3], 1, 3 )
      tmp_tag[j] <- gsub( x = paste( p1[j], "_", p2[j],
                                         tmp_p3, "_", p4[j], sep = "" ),
                              pattern = "\\.", replacement = "" )
      out_mat[j,5] <- tmp_tag[j]
    }
  }
  # Output the newly formatted sequence alignment and the table
  rownames( tmp_ff ) <- out_mat[,5]
  if( ! dir.exists( "../out_RData" ) ) {
    dir.create( "../out_RData")
  }
  write.table( out_mat, file = paste( "../out_RData/", fasta_fnames[i],
                                      ".tsv", sep = "" ),
               sep = "\t", row.names = FALSE, col.names = TRUE, quote = FALSE )
  ape::write.FASTA( x = tmp_ff,
                    file = paste( "../00_raw_data/alignment/", fasta_fnames[i],
                                  "_newIDs.fasta", sep = "" ) )
  # Fix tree file now!
  tmp_fname_check <- gsub( x = fasta_fnames[i], pattern = "aln_",
                           replacement = "" )
  tmp_tname_check <- gsub( x = tree_path, pattern = "..*tree_|_calibnames..*",
                           replacement = "" )
  tmp_tname       <- gsub( x = tree_path, pattern = "..*/|\\.tree",
                           replacement = "" )
  pos_tt          <- which( tmp_tname_check == tmp_fname_check )
  tmp_tt <- ape::read.nexus( file = tree_path[pos_tt] )
  count <- 0
  for( k in out_mat[,6] ){
    count <- count + 1
    pos_match <- which( gsub( pattern = "'", replacement = "",
                              x = tmp_tt$tip.label ) == k )
    tmp_tt$tip.label[pos_match] <- out_mat[count,5]
  }
  # Output fixed tree!
  ape::write.tree( phy = tmp_tt, file = paste( "../00_raw_data/trees/",
                                               tmp_tname[pos_tt],
                                               "_newIDs.tree", sep = "" ) )
}



