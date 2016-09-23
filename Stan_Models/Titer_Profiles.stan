data {
    int<lower=0> N;                 		 ## Number of Data Points
    int<lower=0> N_IE;              		 ## Number of Infection Experiments 
    int<lower=0> N_VL;              		 ## Number of Virus Lineages
    int<lower=0> N_CIT;              		 ## Number of Citations
    int<lower=0> N_BS;              		 ## Number of Bird Species
    int<lower=0> N_VL_BS;            		 ## Number of Virus Lineage:Bird Species
    real Titer[N];                  		 ## Titer, with a lower bound of 0 
    real censoring[N]; 				 ## Censored Titer                              
    real<lower=0> Day[N];            		 ## Day
    real<lower=0> Daysq[N];            		 ## Quadratic Day
    real<lower=0> LD[N];			 ## Log Dose
    int<lower=1, upper=N_IE> IE[N];   		 ## Infection Experiment
    int<lower=1, upper=N_VL> VL[N];   		 ## Virus Lineage (Here meaning Genotype: NY99 vs WN02)
    int<lower=1, upper=N_CIT> CIT[N]; 		 ## Citation
    int<lower=1, upper=N_BS> BS[N];  		 ## Bird Species
    int<lower=1, upper=N_VL_BS> VL_BS[N];  	 ## Virus Lineage : Bird Species
    vector[N] SS;				 ## Sample Size 

}

parameters {
  ### Fixed Effects: Intercepts
    vector[N_VL] alpha;   			 ## Virus Lineage                
    	           
  ### Fixed Effects: Slopes
    vector[N_VL] beta_Day1;                 	 ## Day 1 Virus Lineage
    vector[N_VL] beta_Day2;                	 ## Quadratic Day Virus Lineage
    real beta_LD;                                ## Log Dose

  ### Random Effects
    real alpha_IE_r[N_IE];  			 ## Infection Experiment (Intercept)             
    real alpha_CIT_r[N_CIT];   	                 ## Citation (Intercept)
    real alpha_VLBS_r[N_VL_BS];                  ## Virus Lineage : Bird Species (Intercept)
    real beta_Day1_VLBS_r[N_VL_BS];              ## Virus Lineage : Bird Species (Day)
    real beta_Day2_VLBS_r[N_VL_BS];		 ## Virus Lineage : Bird Species (Quadratic Day)

  ### Random Effect Variances
    real<lower=0> sigmasq_alpha_IE;   	 	 ## Infection Experiment (Intercept)
    real<lower=0> sigmasq_alpha_CIT;             ## Citation (Intercept)
    real<lower=0> sigmasq_alpha_VLBS;            ## Virus Lineage : Bird Species (Intercept)
    real<lower=0> sigmasq_beta_Day1_VLBS;  	 ## Virus Lineage : Bird Species (Day)
    real<lower=0> sigmasq_beta_Day2_VLBS;        ## Virus Lineage : Bird Species (Quadratic Day)

  ### gamma parameter
    real<lower=0> sigma;     
}

transformed parameters {
    real<lower=0> sigma_alpha_IE;             
    real<lower=0> sigma_alpha_CIT;               
    real<lower=0> sigma_alpha_VLBS;              
    real<lower=0> sigma_beta_Day1_VLBS;
    real<lower=0> sigma_beta_Day2_VLBS;
  
    vector[N] mu;                                ## Linear predictor for titer
    
  ### Square root
    sigma_alpha_IE = sqrt(sigmasq_alpha_IE);
    sigma_alpha_CIT = sqrt(sigmasq_alpha_CIT);
    sigma_alpha_VLBS = sqrt(sigmasq_alpha_VLBS);
    sigma_beta_Day1_VLBS = sqrt(sigmasq_beta_Day1_VLBS);
    sigma_beta_Day2_VLBS = sqrt(sigmasq_beta_Day2_VLBS);
    
  for (i in 1:N) {
    mu[i] = 
     alpha[VL[i]] +       			 ## intercept for each virus lineage
     alpha_IE_r[IE[i]] +                         ## random intercept for each infection experiment
     alpha_CIT_r[CIT[i]] +                       ## random intercept for each citation
     alpha_VLBS_r[VL_BS[i]] +                    ## random intercept for each bird species 
     (beta_Day1[VL[i]] +                         ## day slope for each virus lineage
     beta_Day1_VLBS_r[VL_BS[i]]) * 	 	 ## random component of day slope for each VL BS                       
     Day[i] +                                    ## across days
     (beta_Day2[VL[i]] +                         ## quadratic day slope for each virus lineage
     beta_Day2_VLBS_r[VL_BS[i]]) *		 ## quadratic day slope for each bird species  	       
     Daysq[i] +                                  ## across days squared
     beta_LD * LD[i];    			 ## log dose
   }
      
}

model {

 ### diffuse priors for each parameter (fixed effects)  
   alpha ~ cauchy(0, 2.5);			 
   beta_Day1 ~ cauchy(0, 5);
   beta_Day2 ~ cauchy(0, 5);
   beta_LD ~ cauchy(0, 5);					
 
 ### variances
   sigmasq_alpha_IE ~ inv_gamma(1.0E-2, 1.0E-2);        	 
   sigmasq_alpha_CIT ~ inv_gamma(1.0E-2, 1.0E-2);	  	 
   sigmasq_alpha_VLBS ~ inv_gamma(1.0E-2, 1.0E-2);             
   sigmasq_beta_Day1_VLBS ~ inv_gamma(1.0E-2, 1.0E-2);
   sigmasq_beta_Day2_VLBS ~ inv_gamma(1.0E-2, 1.0E-2);
   sigma ~ inv_gamma(1.0E-3, 1.0E-3);

  for (l in 1:N_IE) {
    alpha_IE_r[l] ~ normal(0, sigma_alpha_IE);      
  }

  for (n in 1:N_CIT) {
    alpha_CIT_r[n] ~ normal(0, sigma_alpha_CIT);
  }

  for (q in 1:N_VL_BS) {
    alpha_VLBS_r[q] ~ normal(0, sigma_alpha_VLBS);      
    beta_Day1_VLBS_r[q] ~ normal(0, sigma_beta_Day1_VLBS);
    beta_Day2_VLBS_r[q] ~ normal(0, sigma_beta_Day2_VLBS);
  }

 ### Censored Model
  for (i in 1:N) {
      	if (censoring[i] == 0) Titer[i] ~ lognormal(mu[i], sigma/SS[i]);
    else {
         target += lognormal_lcdf(Titer[i] | mu[i], sigma/SS[i]);
         }
  }
   
}