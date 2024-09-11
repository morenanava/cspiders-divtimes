#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------#
# SET WORKING DIRECTORY #
#-----------------------#
library( rstudioapi ) 
# Get the path to current open R script and find main dir "00_Gene_filtering"
path_to_file <- getActiveDocumentContext()$path
script_wd <- paste( dirname( path_to_file ), "/", sep = "" )
wd <- gsub( pattern = "/scripts", replacement = "", x = script_wd )
setwd( wd )

#--------------#
# LOAD OBJECTS #
#--------------#
# Load tree
raw_tt_paths <- list.files( path = "../../00_data_formatting/00_raw_data/trees/",
                            pattern = "newIDs.tree", full.names = TRUE )
raw_tt <- vector( mode = "list", length = length( raw_tt_paths ) )
names( raw_tt ) <- gsub( x = raw_tt_paths,
                         pattern = "..*/|_calibnames_newIDs.tree", 
                         replacement = "" )
for( i in 1:length(raw_tt) ){
  raw_tt[[ i ]] <- ape::read.tree( file = raw_tt_paths[i] )
}

#-------#
# TASKS #
#-------#
# 1. Get an estimate of the calibration set for the root to have the 
#    time for the speciation event at the root, what we will use 
#    to estimate the mean evolutionary rate later. We have a soft-bound
#    calibration (minimum age = 514 |maximum age = 636.1 Ma), and so we
#    will use the mean value as an approximated mean root age in 100Mya
root_age <- mean( c( 5.14, 6.361) ) # 5.7505 (x100 Mya)

# 2. Find tree height. You can use the function `phytools::nodeHeights` to
#    calculateall the heights of the tree. Then, we can extract the maximum
#    height calculated, which will correspond to the length from the root to
#    the highest tip.
tree_height <- vector( mode = "numeric", length( raw_tt ) )
names( tree_height ) <- names( raw_tt )
for( i in 1:length(raw_tt) ){
  tree_height[i] <- max( phytools::nodeHeights( raw_tt[[ i ]] ) )
}
# tree_noeurycyde     tree_noUCEs     tree_supaln 
# 1.590023        1.266268        1.590023

# 3. Estimate mean rate based on the two different time units
#
#    tree_height = mean_rate x root_age --> mean_rate = tree_height / root_age
mean_rate <- tree_height / root_age
## subst/site per time unit
# tree_noeurycyde     tree_noUCEs     tree_supaln 
# 0.2765017       0.2202014       0.2765017

# If we want to know the mean rate in subst/site/year, we apply the time unit. We
# We should get the same estimate regardless the time unit used:
#
# Time unit 100 May (10^8y): 0.2765017 subst/site/1e+08 = 2.77e-10 subst/site/year

# 4. Now, we can build the gamma distribution given that we now have an estimate
#    for the mean rate. We will also use `alpha = 2` as we will start with a 
#    vague distribution. Nevertheless, if you were very sure about the mean 
#    rate, you could build a more constraint prior.
#
#    mean_Gamma = mean_rate = alpha / beta --> beta = alpha / mean_rate
alpha <- 2
beta  <- alpha/mean_rate
# tree_noeurycyde     tree_noUCEs             tree_supaln 
#  7.233229~7.2        9.082595~9.1        7.233229~7.2

# We can plot this distribution
par( mfrow = c( 1, 3 ) )
curve( dgamma( x, shape = 2, rate = beta[1] ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,7.2) " ), 
        col = "black", lty = 1, box.lty = 2 )
curve( dgamma( x, shape = 2, rate = beta[2] ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,9.1) " ), 
        col = "black", lty = 1, box.lty = 2 )
curve( dgamma( x, shape = 2, rate = beta[3] ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,7.2) " ), 
        col = "black", lty = 1, box.lty = 2 )

# 5. Plot the gamma distributions
if ( ! dir.exists( "out_RData" ) ){
  dir.create( "out_RData" )
}
pdf( file = "out_RData/gamma_dists.pdf", paper = "a4" )
par( mfrow = c( 3, 1 ) )
curve( dgamma( x, shape = 2, rate = beta[1] ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,7.2) " ), 
        col = "black", lty = 1, box.lty = 2 )
title( main = "noeurycyde")
curve( dgamma( x, shape = 2, rate = beta[2] ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,9.1) " ), 
        col = "black", lty = 1, box.lty = 2 )
title( main = "noUCEs" )
curve( dgamma( x, shape = 2, rate = beta[3] ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,7.2) " ), 
        col = "black", lty = 1, box.lty = 2 )
title( main = "supaln")
dev.off()

