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

#-------#
# TASKS #
#-------#
# 0. Load tree
tt <- ape::read.tree( "../../00_data_formatting/00_raw_data/trees/tree_bl.tree" )
  
# 1. Find tree height. You can use the function `phytools::nodeHeights` to  
#    calculate all the heights of the tree. Then, we can extract the maximum 
#    height calculated, which will correspond to the length from the root to
#    the highest tip.
theight <- max( phytools::nodeHeights( tt ) )  # 2.900885

# 2. Get an estimate for the root age based on the age constraint. 
#    We have a soft-bound calibration (minimum age = 3,225 Ma 
#    |maximum age = 4,520 Ma), and so we will use the mean value as an
#    approximated mean root age in (time unit = 1,000Ma = 1 Ga ).
root_age <- mean( c( 3.347, 4.520 ) ) # 3.9335 Ga

# 3. Estimate mean evolutionary rate based:
#    tree_height = mean_rate x root_age --> mean_rate = tree_height / root_age
mean_rate <- theight / root_age # 0.7374819

# If we want to know the mean rate in subst/site/year, we just need to rescale
# the time unit (years instead of Ga):
#
# Time unit 1Ga (1e+09y): 0.7374819 subst/site/1e+09 = 7.37e-10 subst/site/year

# 4. Now, we can build the gamma distribution given that we now have an estimate
#    for the mean rate. We will use `alpha = 2` to account for the uncertainty
#    about the root age. Nevertheless, if you were very sure about the mean 
#    rate, you could build a narrower prior.
#
#    mean_Gamma = mean_rate = alpha / beta --> beta = alpha / mean_rate
alpha <- 2
beta  <- alpha/mean_rate # 2.711931
# We can approximate our beta value to (~2.7)
beta <- 2.7

# Next, let's plot the gamma distribution that we shall use for the rate prior
curve( dgamma( x, shape = 2, rate = beta ), from = 0, to = 3, col = "black" )
legend( "topright", legend = c( "G(2,2.7)" ), 
        col = "black", lty = 1, box.lty = 2 )

# Save a file with the rate prior
if ( ! dir.exists( "out_RData" ) ){
  dir.create( "out_RData" )
}
pdf( file = "out_RData/gamma_dists.pdf", paper = "a4" )
curve( dgamma( x, alpha, rate = beta ), from = 0, to = 2, col = "black" )
legend( "topright", legend = c( "G(2,2.7)" ), 
        col = "black", lty = 1, box.lty = 2 )
dev.off()

