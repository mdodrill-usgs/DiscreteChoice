###############################################################################
#                                                                       Feb '19
#          Discrete Choice Stan Model - RBT Prey Selection
#
#  Notes:
#  * Changed the order of the prey to; c(2,4,5,7,19,21,1)
#  * Dropped the first size bin for all prey 
#    - changes to formatter, and subsequent code
#  * 
#  * Added multi-level R2 to the model
#  * Latest version has no worms
#
###############################################################################
rm(list = ls(all = TRUE))
library(dplyr)
library(rstan)
library(here)
# here()

rstan_options(auto_write = TRUE)
options(mc.cores = parallel::detectCores())
# Sys.setenv(LOCAL_CPPFLAGS = '-march=native')  # can't use this with the larger model...Crash!

#-----------------------------------------------------------------------------#

source(paste0(getwd(),"/Code/Functions.R"))

#-----------------------------------------------------------------------------#
# length
# data.in = model.set.up(model.name = "Length")
data.in = model.set.up.no.worms(model.name = "Length")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", 
            "spsz_eta", "spsz_st_eta", "spsz_ind_eta",
            "beta_f_sz_int")


ni = 2000
nt = 1
nb = 1000
nc = 3

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)



save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_Length_2000iter_with_All_RE.RData")

#-----------------------------------------------------------------------------#
# Width
# data.in = model.set.up(model.name = "Width")
data.in = model.set.up.no.worms(model.name = "Width")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", 
            "spsz_eta", "spsz_st_eta", "spsz_ind_eta",
            "beta_f_sz_int")


ni = 2000
nt = 1
nb = 1000
nc = 3

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)

save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_Width_2000iter_with_All_RE.RData")


#-----------------------------------------------------------------------------#
# Area
# data.in = model.set.up(model.name = "Area")
data.in = model.set.up.no.worms(model.name = "Area")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", 
            "spsz_eta", "spsz_st_eta", "spsz_ind_eta",
            "beta_f_sz_int")

ni = 2000
nt = 1
nb = 1000
nc = 3

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)

save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_Area_2000iter_with_All_RE.RData")


#-----------------------------------------------------------------------------#
# Mass
# data.in = model.set.up(model.name = "Mass")
data.in = model.set.up.no.worms(model.name = "Mass")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", 
            "spsz_eta", "spsz_st_eta", "spsz_ind_eta",
            "beta_f_sz_int")


ni = 2000
nt = 1
nb = 1000
nc = 3

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)

save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_Mass_2000iter_with_All_RE.RData")


#-----------------------------------------------------------------------------#
# set up a parallel section....? Not yet




#-----------------------------------------------------------------------------#
# End