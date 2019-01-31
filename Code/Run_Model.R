###############################################################################
#                                                               Jan 19 :( :( :(
#          Discrete Choice Stan Model - RBT Prey Selection
#
#  Notes:
#  * Changed the order of the prey to; c(2,4,5,7,19,21,1)
#  * Dropped the first size bin for all prey 
#    - changes to formatter, and subsequent code
#  * Moved all of the files to new folder "STAN_Mod_7"
#  * Added multi-level R2 to the model
#
###############################################################################
# setwd('C:/Users/mdodrill/Desktop/FB_DOWN/Analysis/DRIFT_DOWNSTREAM')
setwd('U:/Remote')
rm(list = ls(all = TRUE))
library(dplyr)
library(reshape2)

#-----------------------------------------------------------------------------#
# read in the drift data

tmp.A = read.table(paste0(getwd(), "/Data/Drift_A.txt"), sep = "\t", header = TRUE)
colnames(tmp.A)[3:80] = as.character(1:78) 

tmp.A.w = read.table(paste0(getwd(), "/Data/Drift_A_w.txt"), sep = "\t", header = TRUE)


