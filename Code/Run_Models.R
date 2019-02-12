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
data.in = model.set.up(model.name = "Length")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", #"spsz_eta", "spsz_ind_eta",
            "beta_f_sz_int")


# ni = 2000
ni = 10
nt = 1
# nb = 1000
nb = 5
# nc = 3
nc = 1

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)



# save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_ADDNAME.RData")

#-----------------------------------------------------------------------------#
# Width
data.in = model.set.up(model.name = "Width")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", #"spsz_eta", "spsz_ind_eta",
            "beta_f_sz_int")


# ni = 2000
ni = 10
nt = 1
# nb = 1000
nb = 5
# nc = 3
nc = 1

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)

# save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_ADDNAME.RData")


#-----------------------------------------------------------------------------#
# Area
data.in = model.set.up(model.name = "Area")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", #"spsz_eta", "spsz_ind_eta",
            "beta_f_sz_int")


# ni = 2000
ni = 10
nt = 1
# nb = 1000
nb = 5
# nc = 3
nc = 1

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)

# save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_ADDNAME.RData")


#-----------------------------------------------------------------------------#
# Mass
data.in = model.set.up(model.name = "Mass")

list2env(data.in, globalenv())

params <- c("beta_sz", "mu_sp", "sig_sp",  # the paramaters
            "mu_sp_all", "sig_sp_all", # the transformed paramaters
            "sig_mu", "mu_lprod", "sig_lprod",
            "gamma",
            "n_eprod_a",# "n_tmp_sum",
            # "p", "p_a",
            "vRE", "vRE_sp", "vRE_sp_fsz", "vRE_sp_psz_fsz", # R^2 parms
            "mu_beta_sz", "sig_spsz", "sig_spsz_ind", "sig_spsz_st", #"spsz_eta", "spsz_ind_eta",
            "beta_f_sz_int")


# ni = 2000
ni = 10
nt = 1
# nb = 1000
nb = 5
# nc = 3
nc = 1

Sys.time()

m <- stan_model(paste0(getwd(), "/Stan_Code/",'Model_Ind_v1_with_R2.stan'))
fit <- sampling(m, data = data.in, 
                pars = params,
                control = list(max_treedepth = 12),
                # control = list(max_treedepth = 14, adapt_delta = .925),
                chains = nc, thin = nt, iter = ni, warmup = nb)

# save.image("U:/Desktop/Fish_Git/DiscreteChoice/working_runs/Model_Ind_v1_with_R2_ADDNAME.RData")


#-----------------------------------------------------------------------------#
# set up a parallel section....? Not yet




#-----------------------------------------------------------------------------#
# End