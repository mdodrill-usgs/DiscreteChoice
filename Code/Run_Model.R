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
rm(list = ls(all = TRUE))
library(dplyr)
# library(reshape2)
library(rstan)

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
Sys.setenv(LOCAL_CPPFLAGS = '-march=native')

#-----------------------------------------------------------------------------#
# read in the drift data

tmp.A = read.table(paste0(getwd(), "/Data/Drift_A.txt"), sep = "\t", header = TRUE)
colnames(tmp.A)[3:80] = as.character(1:78) 

# only by taxa, not size
tmp.A.w = read.table(paste0(getwd(), "/Data/Drift_A_w.txt"), sep = "\t", header = TRUE)

identical(tmp.A[,1:2], tmp.A.w[,1:2]) # these should match
#-----------------------------------------------------------------------------#
# read in the diet data

y.tmp2 = read.table(paste0(getwd(), "/Data/Diet.txt"), sep = "\t", header = TRUE)
colnames(y.tmp2)[6:83] = as.character(1:78) 

# only by taxa, not size
w.tmp2 = read.table(paste0(getwd(), "/Data/Diet_w.txt"), sep = "\t", header = TRUE)

identical(y.tmp2[,1:3], w.tmp2[,1:3]) # these should match
#-----------------------------------------------------------------------------#
# format/match/check data

drift.ts = paste(tmp.A[,1], tmp.A[,2])
diet.ts = sort(unique(paste(y.tmp2[,1], y.tmp2[,2])))

# subset only the drift that matches the diet
A = as.matrix(tmp.A[match(diet.ts, drift.ts),3:ncol(tmp.A)])
w.a = as.matrix(tmp.A.w[match(diet.ts, drift.ts),3:ncol(tmp.A.w)])
Nst = as.numeric(nrow(A))  # number of site & trips

# diet
y.tmp = y.tmp2[order(paste(y.tmp2[,1], y.tmp2[,2])),]
w.tmp = w.tmp2[order(paste(w.tmp2[,1], w.tmp2[,2])),]

identical(y.tmp[,1:3], w.tmp[,1:3])
y.in = as.matrix(y.tmp[,6:ncol(y.tmp)])
w.in = as.matrix(w.tmp[,4:ncol(w.tmp)])

#-----------------------------------------------------------------------------#
# format some variables for the model
Nind = nrow(y.in)   # Number of individuals
Nsp = 7             # number of prey taxa 
upper = c(7, 15, 10, 13, 8, 12, 20) - 1
Nspsz = sum(upper)  # number of taxa and sizes

spp = sort(as.numeric(rep(c(1:Nsp), upper)))
# new varible for the fixed paramaters in the drift portion of the model 
spp2 = sort(as.numeric(rep(c(1:Nsp), upper-1)))

# fix all this shit...........
idx = c(rep(1, upper[1]),
        rep(upper[1] + 1, upper[2]),
        rep(upper[1] + upper[2] + 1, upper[3]),
        rep(upper[1] + upper[2] + upper[3] + 1, upper[4]),
        rep(upper[1] + upper[2] + upper[3] + upper[4] + 1, upper[5]),
        rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + 1, upper[6]),
        rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6] + 1, upper[7]))

idx2 = c(rep(upper[1], upper[1]),
         rep(upper[1] + upper[2], upper[2]),
         rep(upper[1] + upper[2] + upper[3], upper[3]),
         rep(upper[1] + upper[2] + upper[3] + upper[4], upper[4]),
         rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5], upper[5]),
         rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6], upper[6]),
         rep(upper[1] + upper[2] + upper[3] + upper[4] + upper[5] + upper[6] + upper[7], upper[7]))

# better way to do this? 
out = NULL
for(i in 1:Nsp){
  out[[i]] = seq(2,upper[i] + 1)  
}

sz = do.call(c,out)

u.idx = unique(idx)
u.idx2 = unique(idx2)

# need these 
Ns = 5
Nt = length(unique(y.tmp[,1]))
trip = as.numeric(y.tmp[,1])
site = as.numeric(y.tmp[,2])

#-----------------------------------------------------------------------------#










