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

#----------------------------#
# DEFINE GLOBAL VARS BY USER #
#----------------------------#
# Name of output calibrated tree file ready to be used by `MCMCtree`.
# Note that the file name will have the following format
# "<your_selected_out_name>_calib_MCMCtree.tree".
#
# NOTE: Make always sure that there is at least one blank line at the 
# end of the this text file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files. This file needs to be
# already in PHYLIP format. Please follow the same format as used in the 
# example tree file provided.
out_name <- c( "anoxphoto" )

# Path to your input text file that allows you to match the flags you have 
# used to label the nodes that are to be calibrated with the calibration you 
# want to use in `MCMCtree` format. If in the same directory,
# you only need to write the name. If in a different directory, please
# type the absolute or relative path to this file.
# The format of this file meets the following criteria:
#
#   - Header.
#   - One row per calibration.
#   - No spaces at all, semi-colon separated,
#   - There are 8 columns:
#       - Name of calibration (no spaces!).
#       - Name of tip 1 that leads to MRCA (no spaces!).
#       - Name of tip 2 that leads to MRCA (no spaces!).
#       - Calibration in `MCMCtree` format. Please check the PAML
#         documentation for more details on this!
# 
# E.g. 1: row in this text file to calibrate node "LUCA"
#
# ```
# # With column 8 only
# LUCA;GCA_000008085;GCA_000021645;'B(3.347,4.520,1e-300,1e-300)'
# ```
#
# If in doubt, please follow the same format as used in the calibrations file 
# used for this analysis
path_textconv_ArchExclDPANN <- c( "../00_raw_data/calibs/Calibs_anoxphoto_withArchExclDPANN.txt" )
calibrations_ArchExclDPANN <- read.table( file = path_textconv_ArchExclDPANN,
                                          stringsAsFactors = FALSE, sep = ";",
                                          blank.lines.skip = TRUE,
                                          header = TRUE,
                                          colClasses = rep( "character", 4 ) )
path_textconv_noArchExclDPANN <- c( "../00_raw_data/calibs/Calibs_anoxphoto.txt" )
calibrations_noArchExclDPANN <- read.table( file = path_textconv_noArchExclDPANN,
                                            stringsAsFactors = FALSE, sep = ";",
                                            blank.lines.skip = TRUE,
                                            header = TRUE,
                                            colClasses = rep( "character", 4 ) )
calibrations <- keep_indexes <- ind_dup <- nodes_dup <- tt_all <-
  vector( mode = "list", 2 )
calibrations[[ 1 ]] <- calibrations_ArchExclDPANN
calibrations[[ 2 ]] <- calibrations_noArchExclDPANN
names( calibrations ) <- names( keep_indexes )  <- names( ind_dup ) <- 
  names( nodes_dup ) <- names( tt_all ) <- c( "withArchExclDPANN",
                                              "withoutArchExclDPANN" )

# Path to tree
path_tree <- c( "../00_raw_data/trees/tree_nobl.tree" )
for( c in 1:length( calibrations ) ){
  cat( "\n[[ ANALYSING CALIBRATION FILE ", names(calibrations)[c], " ]]\n" )
  #tt_ape <- tt_ape_cals <- ape::read.tree( file = path_tree )
  tt_ape <- ape::read.tree( file = path_tree )
  keep_indexes[[c]] <- matrix( 0, nrow = length(rownames(calibrations[[c]])),
                               ncol = 3 )
  # Generate empty vector with as many entries as nodes in the tree
  tt_ape$node.label <- rep( NA, tt_ape$Nnode )
  for( i in 1:length(rownames(calibrations[[c]])) ){
    ## Build MCMCtree calib
    node_lab <- calibrations[[c]][i,4]
    # Get MRCA for these two tips
    mrca <- ape::getMRCA( phy = tt_ape, tip = c( calibrations[[c]][i,2],
                                                 calibrations[[c]][i,3]) )
    keep_indexes[[c]][i,1] <- mrca-ape::Ntip(tt_ape)
    keep_indexes[[c]][i,2] <- calibrations[[c]][i,1]
    keep_indexes[[c]][i,3] <- paste( calibrations[[c]][i,2], "-",
                                calibrations[[c]][i,3], "-",
                                node_lab, sep = "" )
    
    print(mrca-ape::Ntip(tt_ape))
    # Replace node label accordingly
    tt_ape$node.label[mrca-ape::Ntip(tt_ape)] <- paste0( "[",
                                                         calibrations[[c]][i,1],
                                                         "]", collapse = "" )
  }
  ## Find duplicates
  ind_dup[[c]]   <- which( duplicated(keep_indexes[[c]][,1]) == TRUE )
  nodes_dup[[c]] <- which( keep_indexes[[c]][,1] %in% as.numeric( keep_indexes[[c]][ind_dup[[c]],1] ) )
  keep_indexes[[c]][nodes_dup[[c]],]
  # Save tree with node labels
  tt_all[[c]] <- tt_ape
  
} 

##>> ---
## CHECK for duplicates
ind_dup
## NO duplicates, success!
##>> ---

# Remove "NA" from the labs
for( c in 1:length( calibrations ) ){
  ind_na_bools <- is.na( x = tt_all[[c]]$node.label )
  ind_na       <- which( ind_na_bools == TRUE )
  tt_all[[c]]$node.label[ind_na] <- ""
  # Write PHYLIP header, then the calibrated tree
  writeLines( text = paste( length(tt_all[[c]]$tip.label), " 1", sep = "" ), 
              con = paste( "../00_raw_data/trees/cals_only_",
              names( calibrations )[c], ".tree", sep = "" ) )
  ape::write.tree( phy = tt_all[[c]],
                   file = paste( "../00_raw_data/trees/cals_only_",
                   names( calibrations )[c], ".tree", sep = "" ),
                   append = TRUE )
}

#>> TEST
tt_all[[c]]$node.label[c(1,677,687)]
# plot.phylo(tt_all[[c]])
# nodelabels(text=1:tt_all[[c]]$Nnode,node=1:tt_all[[c]]$Nnode+Ntip(tt_all[[c]]),
#            cex = 0.7, frame = "none" )
#>> SUCCESS!
# PLAN: Load the tree later to then replace with the corresponding
# calibrations following my old script

#---------------------------------#
# READ TREE AND CALIBRATIONS FILE #
#---------------------------------#
# Read treea and get phylip header
# NOTE: Make always sure that there is at least one blank line at the 
# end of the tree file! Otherwise, you will get an error telling you that 
# there is an incomplete final line in these files.
tt_name <- c( "../00_raw_data/trees/cals_only_withArchExclDPANN.tree",
              "../00_raw_data/trees/cals_only_withoutArchExclDPANN.tree" )
for( t in 1:length( tt_name ) ){
  
  #-------------------#
  # Get PHYLIP header #
  #-------------------#
  tt            <- readLines( tt_name[t] )
  phylip.header <- tt[1]
  tt            <- tt2 <- tt3 <- tt[2]
  
  #--------------------------------#
  # REPLACE TAGS WITH CALIBRATIONS #
  #--------------------------------#
  # Replace calibration names with corresponding calibration
  for( j in 1:length(rownames(calibrations[[t]])) ){
    # Get MCMCtree calib that will be replaced
    node_lab <- calibrations[[t]][j,4]
    # Now, `node_lab` will have the calibration in `MCMCtree` format already
    # ready to replace the corresponding flag generated in the previous step
    # for each calibrated node
    #
    # Conditional is used so that the single quotation marks are only kept 
    # in the upper-bound calibration for the root. Inequality calibrations
    # do not require single quotation marks
    tmp_calib <- gsub( x = node_lab, pattern = "\\(..*",
                       replacement = "" )
    tmp_calib <- gsub( x = tmp_calib, pattern = "[0-9]..*",
                       replacement = "" )
    if( tmp_calib == 'B' || tmp_calib == 'U' || tmp_calib == 'L' ){
      tt <- gsub( pattern = paste0("\\[",calibrations[[t]][j,1],"\\]"),
                  x = tt,
                  replacement = paste( "'", node_lab, "'", sep = "" ) )
    }else{ # For cross-braced nodes
      tt <- gsub( pattern = paste0("\\[",calibrations[[t]][j,1],"\\]"),
                  x = tt,
                  replacement = paste( node_lab, sep = "" ) )
    }
    # Copy to visualise in FigTree
    reps <- gsub( x = gsub( x = gsub( x = gsub( x = gsub( x = node_lab,
                                                          pattern = "\\{",
                                                          replacement = "(" ),
                                                pattern = "\\}",
                                                replacement = ")" ), 
                                      pattern = "\\[|\\]", replacement = "" ),
                            pattern = "\\#", replacement = "flag" ),
                  pattern = " ", replacement = "-" )
    # For cross-braced calibrations without fossil
    if( tmp_calib == '#' ){
      reps <- gsub( x = gsub( x = reps, pattern = "\\#", replacement = "flag" ),
                    pattern = "\\]", replacement = "" )
      tt2 <- gsub( pattern = paste0("\\[",calibrations[[t]][j,1],"\\]"),
                   x = tt2,
                   replacement = paste0( "'", reps, "-",
                                         calibrations[[t]][j,1], "'", 
                                         collapse = "" ) )
    }else{ # For the rest of calibrations
      tt2 <- gsub( pattern = paste0("\\[",calibrations[[t]][j,1],"\\]"),
                   x = tt2,
                   replacement = paste0( "'", reps, "-",
                                         calibrations[[t]][j,1], "'",
                                         collapse = "" ) )
    }
    # Generate an uncalibrated tree for CODEML!
    tt3 <- gsub( pattern = paste0("\\[",calibrations[[t]][j,1],"\\]"),
                 x = tt3,
                 replacement = "" )
  }
  
  
  #-------------------------------#
  # WRITE CALIBRATED TREE IN FILE #
  #-------------------------------#
  out_dir <- "../01_inp_data/"
  if( ! dir.exists( "../01_inp_data/" ) ){
    dir.create( "../01_inp_data/" )
  }
  write( x = phylip.header, file = paste( out_dir, out_name, "_",
                                          names( calibrations )[t],
                                          "_calib_MCMCtree.tree", sep = "" ) )
  write( x = tt, file = paste( out_dir, out_name, "_",
                               names( calibrations )[t],
                               "_calib_MCMCtree.tree", sep = "" ),
         append = TRUE )
  write( x = phylip.header, file = paste( "../00_raw_data/trees/",
                                          out_name, "_",
                                          names( calibrations )[t],
                                          "_fordisplay_calib_MCMCtree.tree",
                                          sep = "" ) )
  write( x = tt2, file = paste( "../00_raw_data/trees/", out_name, "_",
                                names( calibrations )[t],
                                "_fordisplay_calib_MCMCtree.tree", sep = "" ),
         append = TRUE )
  # Write an uncalibrated tree only once!
  if( t == 1 ){
    write( x = phylip.header, file = paste( out_dir, out_name,
                                            "_uncalib.tree", sep = "" ) )
    write( x = tt3, file = paste( out_dir, out_name,
                                  "_uncalib.tree", sep = "" ),
           append = TRUE )
  }
  
}

  




