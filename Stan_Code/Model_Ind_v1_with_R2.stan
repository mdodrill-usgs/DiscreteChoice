// Feb 2019
// * The prior version of this was Model_Ind_v1.stan

data {
  int Nspsz;           // Number of taxa & size bins      
  int Nst;             // Number of sites & trips
  int Nsp;             // Number of taxa
  int Nind;            // Number of individuals
  
  int idx_ts_ind[Nind];// index for trip and site of individuals
  int u_idx[Nsp];      // Diet
  int u_idx2[Nsp];     // Diet
  int idx[Nspsz];      // Drift
  int idx2[Nspsz];     // Drift
  
  int idx_first[Nsp];         //Drift index for the first position of each taxa in the size bin taxa matrix (i.e., lprod)  
  int sp_idx2[Nspsz - Nsp];   //Drift new index for each taxa and size bin - 1, for each taxa
  int not_first[Nspsz - Nsp]; //Drift oposite of idx_first
  
  vector[Nspsz] sz;        // sizes - diet
  vector[Nind] fish_sz;    // centered fork length of fish
  real avg_log_len;
  
  int y[Nind, Nspsz];      // diet matrix 
  int w[Nind, Nsp];        // diet matrix - extra
  matrix[Nspsz, Nsp] X;    // dummy
  
  int a[Nst, Nspsz];       // drift matrix
  int w_a[Nst, Nsp];       // drift matrix - extra
  vector[Nsp] emp_a;       // empirical proportions of taxa in the drift (across all samples)
  int sp[Nspsz];
  
  vector<lower = 0> [Nsp]alpha;  // dirichlet params 
} 

parameters {
  // Drift
  real<lower = 0> sig_mu;
  vector<lower = 0>[Nsp] sig_lprod;
  vector[Nspsz - Nsp] mu_lprod;
  matrix[Nst, Nspsz - Nsp]lprod;
  simplex[Nsp] tmp_ps_a[Nst];
  
  // Diet
  vector<lower = 0>[Nsp-1] sig_sp;
  
  real mu_beta_sz;
  real<lower = 0> sz_sig_eta;
  real<lower = 0> mu_sp_sig_eta;
  
  vector[Nsp-1] sp_int;
  
  vector[Nspsz] spsz_eta;
  matrix[Nind, Nspsz]spsz_ind_eta;
  matrix[Nst,Nspsz]spsz_st_eta;
  
  real<lower = 0>sig_spsz;
  real<lower = 0>sig_spsz_ind;
  real<lower = 0>sig_spsz_st; 
  
  real beta_f_sz_int;
  
}

transformed parameters {
  // Drift
  matrix[Nst, Nspsz] eprod_a;
  matrix[Nst, Nspsz] n_eprod_a;
  matrix[Nst, Nspsz] p_a;
  matrix[Nst, Nsp] ps_a;
  matrix[Nst, Nspsz] sum_eprod_a;
  matrix [Nst, Nspsz]fix_lprod;
  
  // Diet
  matrix[Nind, Nspsz] p;
  matrix[Nind, Nsp] ps;
  matrix[Nind, Nspsz] beta1;
  matrix[Nind, Nspsz] eprod_u;
  vector[Nind] sum_eprod_u;
  
  real beta_sz;
  vector[Nsp] mu_sp;
  
  // drift ////////////////////////////////////////////
  for(i in 1:Nst){                 
     ps_a[i,] = tmp_ps_a[i]'; 
  } 
  
  for(i in 1:Nst){
    for(j in 1:Nsp){
      fix_lprod[i,idx_first[j]] = 0.0;
    }
  }
  
  for(i in 1:Nst){
    for(j in 1:Nspsz - Nsp){
      fix_lprod[i,not_first[j]] = lprod[i,j];
    }
  }
  
    for(j in 1:Nspsz){
      eprod_a[,j] = exp(fix_lprod[,j]);
    }
    
    for(i in 1:Nst){       // error if I remove the lines below from the i loop
      for(j in 1:Nspsz){
        sum_eprod_a[i,j] = sum(eprod_a[i,idx[j]:idx2[j]]);
      }
      
      for(j in 1:Nspsz){
        n_eprod_a[i,j] = eprod_a[i,j] / sum_eprod_a[i,j];   
        p_a[i,j] = n_eprod_a[i,j] * ps_a[i,sp[j]];   
      }
    }
    
  // diet ////////////////////////////////////////////
  
  beta_sz = exp(mu_beta_sz);

  for(i in 1:Nsp - 1){
    mu_sp[i] = sp_int[i];
  }
  mu_sp[7] = 0.0;
  
  for(i in 1:Nind){
    for(k in 1:Nspsz){
      beta1[i,k] = dot_product(mu_sp, X[k,]) + spsz_eta[k] + beta_sz * (log(sz[k]) - avg_log_len) +
      spsz_ind_eta[i,k] + spsz_st_eta[idx_ts_ind[i],k] + beta_f_sz_int * fish_sz[i] * (log(sz[k]) - avg_log_len);
    }  
  }
  
  for(k in 1:Nspsz){
    eprod_u[,k] = p_a[idx_ts_ind[],k] .* exp(beta1[,k]);  // where the avail comes in 
  }
     
  for(i in 1:Nind){ 
    sum_eprod_u[i] = sum(eprod_u[i,]);
  }
     
  for(k in 1:Nspsz){
    p[,k] = eprod_u[,k] ./ sum_eprod_u;
  }  
 
  for(i in 1:Nind){    
     for(k in 1:Nsp){
       ps[i,k] = sum(p[i, u_idx[k]:u_idx2[k]]);
     }
  }
}

model {
   // drift ////////////////////////////////////////////
  for(i in 1:Nst){
    tmp_ps_a[i] ~ dirichlet(alpha[]);
  }
  
  sig_mu ~ normal(0, 10);  
  sig_lprod ~ normal(0, 10);
  mu_lprod ~ normal(0, sig_mu); 
    
  for(j in 1:Nspsz - Nsp){
    lprod[,j] ~ normal(mu_lprod[j], sig_lprod[sp_idx2[j]]);
  }
    
  for(i in 1:Nst){
    w_a[i,1:Nsp] ~ multinomial(ps_a[i,]');
    a[i,1:Nspsz] ~ multinomial(p_a[i,]');
  }
  
   // diet ////////////////////////////////////////////
   beta_f_sz_int ~ normal(0,10);
   sig_sp ~ normal(0, 10);
   
   sig_spsz ~ normal(0,10);
   sig_spsz_ind ~ normal(0,10);
   
   sp_int ~ normal(0,10);
   
   spsz_eta ~ normal(0,sig_spsz);
   
   for(i in 1:Nind){
     for(k in 1:Nspsz){
       spsz_ind_eta[i,k] ~ normal(0, sig_spsz_ind);
     }
   }
   
   sig_spsz_st ~ normal(0, 10);
   
   for(i in 1:Nst){
     for(k in 1:Nspsz){
       spsz_st_eta[i,k] ~ normal(0, sig_spsz_st);
     }
   }
   
   sz_sig_eta ~ normal(0, 1);
   mu_beta_sz ~ normal(0, 10);
   
   mu_sp_sig_eta ~ normal(0, 1);
   
  for(i in 1:Nind){    
    y[i,1:Nspsz] ~ multinomial(p[i,]');
    w[i,1:Nsp] ~ multinomial(ps[i,]');
  }
}

generated quantities {
  // full set of parameters
  vector[Nsp] mu_sp_all;
  vector[Nsp] mu_sp_tmp;
  vector[Nsp] sig_sp_all;
  vector[Nsp] sig_sp_tmp;

  // Calculating 'gamma' 
  vector[Nspsz] mu_lprod_all;
  vector[Nspsz] tmp_sum;
  vector[Nspsz] n_tmp_sum;
  vector[Nspsz] new_avail;
  vector[Nspsz] new_tmp;
  vector[Nsp] new_tmp_sum;
  vector[Nsp] gamma;
  
  // R^2
  real vRE_sp_psz_fsz;
  row_vector[Nspsz] tmp_sp_sz_eta;
  matrix[Nind, Nspsz] tmp_vRE;
  real vRE;
  matrix[Nind, Nspsz] tmp_vRE_sp;
  real vRE_sp;
  matrix[Nind, Nspsz] tmp_vRE_sp_fsz;
  real vRE_sp_fsz;
  
  // Calculate the transformed intercepts & their sd
  mu_sp_tmp[Nsp] = 0.0;
  sig_sp_tmp[Nsp] = 0.0;

  for(i in 1:Nsp-1){
    mu_sp_tmp[i] = mu_sp[i];
    sig_sp_tmp[i] = sig_sp[i];
  }

  for(i in 1:Nsp){
    mu_sp_all[i] = mu_sp_tmp[i] - mean(mu_sp_tmp[]);
    sig_sp_all[i] = sig_sp_tmp[i] - mean(sig_sp_tmp[]);
  }
  
  // ----------------------------------- //
  // calculate 'gamma' 
  // fill in the zeros for this
  for(i in 1:Nsp){
    mu_lprod_all[idx_first[i]] = 1;
  }

  // // fill in the values of mu_lprod
  for(i in 1:Nspsz - Nsp){
    mu_lprod_all[not_first[i]] = exp(mu_lprod[i]);
  }

  // // sums for each taxa
  for(i in 1:Nspsz){
    tmp_sum[i] = sum(mu_lprod_all[idx[i]:idx2[i]]);
  }

  // // sums to 1 within taxa
  for(i in 1:Nspsz){
    n_tmp_sum[i] = mu_lprod_all[i] / tmp_sum[i];
  }

  // // scale by empirical proportions obs. in drift.
  // // this vector sums to 1
  for(i in 1:Nspsz){
    new_avail[i] = n_tmp_sum[i] * emp_a[sp[i]];
  }

  // // should this be mu_sp_all? (it doesn't matter, mu_sp_all or mu_sp_tmp) Add the st deviates in?
  for(i in 1:Nspsz){
    new_tmp[i] = new_avail[i] * exp(mu_sp_all[sp[i]] + spsz_eta[i] + beta_sz * (log(sz[i]) - avg_log_len));
  }

  for(i in 1:Nsp){
    new_tmp_sum[i] = sum(new_tmp[u_idx[i]:u_idx2[i]]);
  }

  for(i in 1:Nsp){
    gamma[i] = new_tmp_sum[i] / sum(new_tmp_sum[]);
  }
  
  // ----------------------------------- //
  // multi-level R^2                                  
  vRE_sp_psz_fsz = variance(beta1[]);

  for(i in 1:Nspsz){
    tmp_sp_sz_eta[i] = spsz_eta[i];  //convert this to a row vector
  }

  for(i in 1:Nind){
      tmp_vRE[i,] = tmp_sp_sz_eta[] + spsz_st_eta[idx_ts_ind[i],] + spsz_ind_eta[i,];
  }

  vRE = variance(tmp_vRE[]);

  for(k in 1:Nspsz){
    tmp_vRE_sp[,k] = tmp_vRE[,k] + dot_product(mu_sp, X[k,]);
  }
  
  vRE_sp = variance(tmp_vRE_sp[]);
  
  for(i in 1:Nind){
    for(k in 1:Nspsz){
      tmp_vRE_sp_fsz[i,k] = tmp_vRE[i,k] + dot_product(mu_sp, X[k,]) + beta_f_sz_int * fish_sz[i] * (log(sz[k]) - avg_log_len);
    }
  }

  vRE_sp_fsz = variance(tmp_vRE_sp_fsz[]);
    
}
