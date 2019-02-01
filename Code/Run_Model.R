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
spsz = c(1:Nspsz)

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
# calc the average size of prey in the drift (probably a better way :/)
temp = data.frame(sz = sz,
                  sums = colSums(A))

temp.2 = group_by(temp, sz) %>%
  summarize(my.sum = sum(sums))

tmp.vec = list()

for(i in 1:dim(temp.2)[1]){
  tmp.vec[[i]] = rep(as.character(temp.2[i,1]), times = as.numeric(temp.2[i,2]))
}

all.len = do.call('c', tmp.vec)

avg.log.len = mean(log(as.numeric(all.len)))

#-----------------------------------------------------------------------------#
# need emp_a (drift prop for each taxa across all the data)
tmp.a = data.frame(spp = spp,
                   sums = colSums(A))

tmp.a.2 = group_by(tmp.a, spp) %>%
  summarize(my.sum = sum(sums))

a.tmp.w = colSums(w.a)

tmp.a.all = tmp.a.2$my.sum + a.tmp.w

emp_a_2 = tmp.a.all / sum(tmp.a.all)

# emp.dat2 = data.frame(taxa = name.key$name,
#                       emp = emp_a_2)
#-----------------------------------------------------------------------------#
# create a vector with the same size distribution between bins (the overall size
# dist), but with the same relative diff in abudance between taxa (bunch of code
# for this --> clean up)

# sum the over 8mm into one bin
dat = data.frame(spp = spp,
                 sz = sz,
                 sums = colSums(A))

out = expand.grid(spp = unique(spp), sz = 2:8)
out = out[order(out$spp),]
out$count = 0

big.out = NULL

for(i in 1:Nsp){
  sub = dat[which(dat$spp == i),]
  out.s = out[which(out$spp == i),]
  
  out.s[1:6,3] = sub$sums[1:6]
  
  if(nrow(sub) >= 7){
    out.s[7,3] = sum(sub$sums[7:nrow(sub)])
  }
  big.out[[i]] = out.s
}

out2 = do.call('rbind', big.out)

# get the overall size distribution, across taxa
dist = group_by(out2, sz) %>%
  summarize(my.sum = sum(count))

dist$prop = dist$my.sum / sum(dist$my.sum)

# total count and prop for each taxa
tot.ct = group_by(out2, spp) %>%
  summarize(tot = sum(count))
tot.ct$prop = tot.ct$tot / sum(tot.ct$tot)

out3 = list()

for(i in 1:Nsp){
  out3[[i]] = dist$prop * tot.ct$prop[i]
}

out4 = do.call(c, out3)
sum(out4)

eq_dist = out4
eq_sp = out2$spp
eq_sz = out2$sz
Neqsz = dim(out)[1]

eq_log_len = mean(log(eq_sz))

#-----------------------------------------------------------------------------#
# area calc (see C:\Users\mdodrill\Desktop\FB_DOWN\Analysis\IMAGE\Image_V2.R)
# wrote using 'dput'

fit.dat = structure(list(V1 = c("SIMA", "SIML", "CHIA", "CHIL", "GAM", "NZMS", "LUMB"),
                         V2 = c(0.224460916479209, -1.78662907514783, -1.34743919150295,
                                -2.15070453700794, -1.45525542402561, -0.74863545618833,
                                -1.85570798696119),
                         V3 = c(0.94958327948254, 1.79612553310148, 1.54243803018606,
                                1.67798706868551, 1.89261925347815, 1.65868421232674,
                                1.30725930407769)),
                    row.names = c(NA, -7L), class = "data.frame")

name.key = data.frame(num = c(1:7),
                      abr.name = c("NZMS", "GAM", "SIMA", "SIML",
                                   "CHIA", "CHIL", "LUMB"),
                      name = c("Snails", "Gammarus",
                               "Black Fly Adult", "Black Fly Larva",
                               "Midge Adult", "Midge Larva","Worms"))
tmp.ar2 = list()

fit.c = vector(length = Nsp)

del = vector(length = Nsp)

area.df = list()

for(i in 1:Nsp){
  sub = fit.dat[which(fit.dat$V1 == name.key[i,2]),]
  
  fit.c[i] = sub[,3]
  
  tmp.sz = sz[which(spp == i)]
  
  tmp.ar = sub[,2] + sub[,3] * log(tmp.sz)
  
  del[i] = sub[,2] + sub[,3] * avg.log.len  #log area for average size
  
  tmp.ar2[[i]] = exp(tmp.ar)
  
  area.df[[i]] = data.frame(Taxa = name.key[i,3],
                            Size = tmp.sz,
                            area = exp(tmp.ar))
}

# area = do.call('c', tmp.ar2)
area = do.call('rbind', area.df)

delta = vector(length = Nsp - 1)

for(i in 1:Nsp-1){
  # delta[i] = (fit.c[i] - fit.c[7]) * exp(avg.log.len)       # old
  delta[i] = (del[i] - del[7])
}

#-----------------------------------------------------------------------------#
alpha = rep(.2, Nsp)

X <- as.matrix(model.matrix(~ as.factor(spp) - 1))   

idx_first = unique(idx)  # position of the first bin for each taxa
tmpper = seq(1:Nspsz)
not_first = tmpper[which(!tmpper %in% idx_first)]
#-----------------------------------------------------------------------------#
params <- c("beta_sz", "mu_sp", "sig_sp",# "beta_sp_st",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            # "gamma_st", "gamma_sz", "ratio_st",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            "p", "p_a",
            "vRE_1", "vRE_2",
            "vRE", "vRE_sp", "vRE_sz", "vRE_sp_sz", 
            "mu_beta_sz", "spsz_eta", "spsz_st_eta", "sig_spsz", "sig_spsz_st")

data.in = list(Nspsz = Nspsz, Nst = Nst, Nsp = Nsp, Neqsz = Neqsz, a = A, sp = spp, sp_idx2 = spp2,
               w_a = w.a, idx = idx, idx2 = idx2, idx_first = idx_first, not_first = not_first,
               alpha = alpha, a_Nsz = upper, # end avail  Nind = Nind,
               Ns = Ns, Nt = Nt, trip = trip, site = site,
               # eq_dist = eq_dist, eq_sp = eq_sp, eq_sz = eq_sz,  
               y = y.in, X = X, Nsp = Nsp,
               emp_a = emp_a_2,
               # sz = sz-3.5, u_Nsz = Nsz, w = w.in, u_idx = u.idx, u_idx2 = u.idx2, c = fit.c,  delta = delta, 
               sz = sz, u_Nsz = Nsz, w = w.in, u_idx = u.idx, u_idx2 = u.idx2,
               avg_log_len = avg.log.len, #eq_log_len = eq_log_len,
               # sz = sz, u_Nsz = Nsz, w = w.in, u_idx = u.idx, u_idx2 = u.idx2,
               # area = area, u_Nsz = Nsz, w = w.in, u_idx = u.idx, u_idx2 = u.idx2,
               spsz = spsz)


## MCMC settings
ni <- 2000
nt <- 2
nb <- 500
nc <- 3

fit <- stan("Model_7_V4.stan", 
            data = data.in, 
            pars = params,
            control = list(max_treedepth = 12),
            # control = list(max_treedepth = 14, adapt_delta = .925),
            chains = nc, thin = nt, iter = ni, warmup = nb)


