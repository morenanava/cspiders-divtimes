###############################################################################
# This script, created for Nava et al 2024 creates a range through time plot  #
# of pycnogonid fossils based on data from the PBDB which it loads, adds taxa #
# to, and then plots. If you run the script line by line, it will create a    #
# PDF in the root directory, which was then modified in Inkscape for the      #
# figure used in publication                                                  #
###############################################################################

# Clean environment
rm( list = ls( ) )

# If packages aren't installed, install them, then load them
packages <- c( "palaeoverse", "rstudioapi" )
if( length(packages[!packages %in% installed.packages()[,"Package"]]) > 0 ){
  install.packages( packages[!packages %in% installed.packages()[,"Package"]] )
}
library( "palaeoverse" )
library( "rstudioapi" )

# Get the path to current open R script and find main dir
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
# If you have had issues installing `rstudioapi`, you may want to manually
# set the working directory to the path of the root folder (i.e. the one in
# which this file is saved) setwd("/home/location/of/files/")

# Load R in-house functions
source( "Functions_RTTplot.R" )

###############################################################################
# The code below will load PBDB data, add two taxa, and then create a PDF in  #
# our root directory                                                          #
###############################################################################

# Load PBDB data from file, exported July 24 from the PBDB
#occurenceData<-read.csv(here("pycno_pbdb_data1.csv"), skip = 19)
occurenceData <- read.csv( file = "../data/pycno_pbdb_data1.csv", skip = 19 )
# Filter out taxa that cannot be identified to species level
occurenceDataSpecies <- occurenceData[occurenceData$accepted_rank=="species",]
# For the graphs in this paper, we would like to group these by whether
# they are stem pycnogonids, or crown
occurenceDataSpecies$taxonomy <- NA
# For these taxa that is relatively easy - those that are Mesozoic are crown,
# those that are Palaeozoic are stem - apply correct labels:
gtmin <- which( occurenceDataSpecies$min_ma>252 )
ltmin <- which( occurenceDataSpecies$min_ma<252 )
occurenceDataSpecies[gtmin,]$taxonomy <- "Stem"
occurenceDataSpecies[ltmin,]$taxonomy <- "Pantopoda"

# Add Cambropycnogon, which is not in PBDB, but is of note
Cambropycnogon <- occurenceDataSpecies[1,]
Cambropycnogon <- Cambropycnogon[-1,]
Cambropycnogon[1,] <- NA
Cambropycnogon[1,]$genus <- "Cambropycnogon"
# Based on age of Rehbachiella in PBDB, which is also found at Kinnekull
Cambropycnogon[1,]$max_ma <- 497.0 
# and in line with Agnostus pisiformis Zone in which this was reportedly found
Cambropycnogon[1,]$min_ma <- 496.8 
Cambropycnogon[1,]$accepted_name <- "Cambropycnogon klausmuelleri"
Cambropycnogon[1,]$accepted_rank <- "species"
Cambropycnogon[1,]$ref_pubyr <- 2002
# There is a good deal of debate about whether this is indeed a pycnogonid,
# hence decision to place as incertae sedis
Cambropycnogon[1,]$taxonomy <- "Incertae sedis"
occurenceDataSpecies <- rbind( occurenceDataSpecies, Cambropycnogon )
rm( Cambropycnogon )

# Add Palaeomarachne, also not in PBDB
Palaeomarachne <- occurenceDataSpecies[1,]
Palaeomarachne <- Palaeomarachne[-1,]
Palaeomarachne[1,] <- NA
Palaeomarachne[1,]$genus <- "Palaeomarachne"
# Based on other William Lake fossil (Lunataspis) in PBDB
Palaeomarachne[1,]$max_ma <- 453.0 
Palaeomarachne[1,]$min_ma <- 445.2
Palaeomarachne[1,]$accepted_name <- "Palaeomarachne granulata"
Palaeomarachne[1,]$accepted_rank <- "species"
Palaeomarachne[1,]$ref_pubyr <- 2002
#Also Incertae sedisto reflect previous work by this team
Palaeomarachne[1,]$taxonomy  <- "Incertae sedis"
occurenceDataSpecies <- rbind( occurenceDataSpecies, Palaeomarachne )
rm( Palaeomarachne )

# Create our plot
if( dir.exists( "../plots/" ) == FALSE ){
  dir.create( "../plots/" )
}
#pdf(here("Pycnogonid_fossils.pdf"), width = 8, height = 5.5)
pdf( file = "../plots/Pycnogonid_fossils.pdf", width = 9, height = 5.5 )
mar.default <- c( 5,4,4,2 ) + 0.1
par( mar = mar.default + c( 2, 0, 0, 0 ) )
tax_range_time_RJG( occdf = occurenceDataSpecies, name = "genus", plot = TRUE,
                    by = "FAD", group_by = "taxonomy",
                    title = "The pycnogonid fossil record", plot_names = TRUE,
                    x_pad = 10, x_range = c( 520, 10 ),
                    # Values that make it possible to have the labels at the
                    # top and at the bottom close to the plot axes
                    y_range = c( 1, 12 ), 
                    # if xlabline = 3, x axis is too close to numbers
                    xlabline = 4 )
dev.off( )
