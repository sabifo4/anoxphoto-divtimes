#-------------------#
# CLEAN ENVIRONMENT #
#-------------------#
rm( list = ls( ) )

#-----------------------------------------------#
# LOAD PACKAGES, FUNCTIONS, AND SET ENVIRONMENT #
#-----------------------------------------------#
# This package lets you find automatically the path to a specific location
# in your file structure
# If you have not installed this package, you will need to install it. 
# You can uncomment the following line to do this:
#install.packages( "rstudioapi" )
library( rstudioapi )
path_to_file <- getActiveDocumentContext()$path
wd <- paste( dirname( path_to_file ), "/", sep = "" )
setwd( wd )
# Set the home directory to where this script is saved
home_dir <- gsub( pattern = "scripts/", replacement = "", x = wd )
# Load the file with all the functions used throughout this script
source( file = "../../../src/Functions.R" )

#-------------------------------------------------------------#
# DEFINE GLOBAL VARIABLES -- modify according to your dataset #
#-------------------------------------------------------------#
# First, we will define the global variables that match the settings in our 
# analysis.

# Read trees
## tt1: tree with node labels as of MCMCtree --> there are two
##      `node_tree_[1-2].tree` files but they are the same, so we will use
##      just one!
## tt2: tree with node labels including PAML calibration and node name
tt1 <- ape::read.tree( file = "../sum_analyses/00_prior/node_tree_1.tree" )
tt2 <- ape::read.tree( file = "../../../00_data_formatting/00_raw_data/trees/anoxphoto_withArchExclDPANN_fordisplay_calib_MCMCtree.tree" )
tt3 <- ape::read.tree( file = "../../../00_data_formatting/00_raw_data/trees/anoxphoto_withoutArchExclDPANN_fordisplay_calib_MCMCtree.tree" )
# Create directory to save files with calibration info
if( ! dir.exists( "../calib_files" ) ){
  dir.create( "../calib_files" ) 
}
# Save the headers first
write( x = "Calib;node;Prior",
       file = "../calib_files/Calibnodes_withArchExclDPANN.csv" )
write( x = "Calib;node;Prior",
       file = "../calib_files/Calibnodes_withArchExclDPANN_margVScalib.csv" )
write( x = "Calib;node;Prior",
       file = "../calib_files/Calibnodes_withoutArchExclDPANN.csv" )
write( x = "Calib;node;Prior",
       file = "../calib_files/Calibnodes_withoutArchExclDPANN_margVScalib.csv" )
count_margvscalib <- 0
## TREE 1
# Now, populate the output file with the info for the rest of the nodes
for( i in 1:length(tt2$node.label) ){
  
  if( tt2$node.label[i] != "" ){
    tmp_name <- gsub( x = gsub( x = gsub( x = tt2$node.label[i],
                                          pattern = "..*\\)-",
                                          replacement = "" ), pattern = "'",
                                replacement =  "" ),
                      pattern = "\\/", replacement = "_" )
    tmp_dist <- gsub( x = gsub( x = tt2$node.label[i], pattern = "\\)-..*",
                                replacement = "\\)" ), pattern = "'",
                      replacement = "" )
    is_notdist <- grep( x = tmp_dist, pattern = "flag" )
    if( length( is_notdist ) != 0 ){
      write( x = paste( tmp_name, ";", tt1$node.label[i], ";", tmp_dist,
                        sep = "" ), 
             file = "../calib_files/Calibnodes_withArchExclDPANN_margVScalib.csv",
             append = TRUE )
      count_margvscalib <- count_margvscalib + 1
    }
    cat( paste( tmp_name, ";", tt1$node.label[i], ";", tmp_dist, "\n",
                sep = "" ) )
    write( x = paste( tmp_name, ";", tt1$node.label[i], ";", tmp_dist,
                      sep = "" ), 
           file = "../calib_files/Calibnodes_withArchExclDPANN.csv",
           append = TRUE )
  }
  
}
# Remove second csv if not needed
if( count_margvscalib == 0 ){
  unlink( x = "../calib_files/Calibnodes_withArchExclDPANN_margVScalib.csv" )
}

## TREE 2
# Now, populate the output file with the info for the rest of the 
count_margvscalib <- 0
for( i in 1:length(tt3$node.label) ){
  
  if( tt3$node.label[i] != "" ){
    tmp_name <- gsub( x = gsub( x = gsub( x = tt3$node.label[i],
                                          pattern = "..*\\)-",
                                          replacement = "" ), pattern = "'",
                                replacement =  "" ),
                      pattern = "\\/", replacement = "_" )
    tmp_dist <- gsub( x = gsub( x = tt3$node.label[i], pattern = "\\)-..*",
                                replacement = "\\)" ), pattern = "'",
                      replacement = "" )
    is_notdist <- grep( x = tmp_dist, pattern = "flag" )
    if( length( is_notdist ) != 0 ){
      write( x = paste( tmp_name, ";", tt1$node.label[i], ";", tmp_dist,
                        sep = "" ), 
             file = "../calib_files/Calibnodes_withoutArchExclDPANN_margVScalib.csv",
             append = TRUE )
      count_margvscalib <- count_margvscalib + 1
    }
    cat( paste( tmp_name, ";", tt1$node.label[i], ";", tmp_dist, "\n",
                sep = "" ) )
    write( x = paste( tmp_name, ";", tt1$node.label[i], ";", tmp_dist,
                      sep = "" ), 
           file = "../calib_files/Calibnodes_withoutArchExclDPANN.csv",
           append = TRUE )
  }
  
}
# Remove second csv if not needed
if( count_margvscalib == 0 ){
  unlink( x = "../calib_files/Calibnodes_withoutArchExclDPANN_margVScalib.csv" )
}

