data {
    int<lower=0> N;                 		 ## Number of Data Points 
    int<lower=0> N_IE;              		 ## Number of Infection Experiments
    int<lower=0> N_VL;              		 ## Number of Virus Lineages
    int<lower=0> N_CIT;              		 ## Number of Citations
    int<lower=0> N_BS;              		 ## Number of Bird Species
    int<lower=0> N_VL_BS;            		 ## Number of Virus Lineage:Bird Species
                                    
    real<lower=0> Day[N];            		 ## Day
    real<lower=0> LD[N];			 ## Log Dose
    int<lower=1, upper=N_IE> IE[N];   		 ## Infection Experimnet
    int<lower=1, upper=N_VL> VL[N];   		 ## Virus Lineage (Here meaning Genotype: NY99 vs WN02)
    int<lower=1, upper=N_CIT> CIT[N]; 		 ## Citation
    int<lower=1, upper=N_BS> BS[N];  		 ## Bird Species
    int<lower=1, upper=N_VL_BS> VL_BS[N];  	 ## Virus Lineage : Bird Species 

    real<lower=0, upper=1> Survival[N]; 	 ## Survival
    int<lower=0> Samp_Max;              	 ## Max number of Samples
    int<lower=0> Ali_Max;              	 	 ## Max number of Living Birds
    int<lower=0, upper=Samp_Max> N_Samp[N];   	 ## Number of Samples
    int<lower=0, upper=Ali_Max> N_Ali[N];     	 ## Number of Living Birds
}

parameters {
  ### Fixed Effects: Intercepts
    vector[N_VL] alpha;             		 ## Virus Lineage         
             
  ### Fixed Effects: Slopes
    vector[N_VL] beta_Day;                 	 ## Day Virus Lineage
    vector[N_VL] beta_LD;                	 ## Log Dose             
 
  ### Random Effects
    real alpha_IE_r[N_IE];                  	 ## Infection Experiment
    real beta_Day_IE_r[N_IE];                      
    real alpha_CIT_r[N_CIT];  		         ## Citation
    real beta_Day_CIT_r[N_CIT];
    real alpha_VLBS_r[N_VL_BS];                  ## Virus Lineage : Bird Species
    real beta_Day_VLBS_r[N_VL_BS] ;             

  ### Random Effect Varaiances
    real<lower=0> sigmasq_alpha_IE;       	 ## Unique Line
    real<lower=0> sigmasq_beta_Day_IE;    
    real<lower=0> sigmasq_alpha_CIT;             ## Citation
    real<lower=0> sigmasq_beta_Day_CIT;
    real<lower=0> sigmasq_alpha_VLBS;            ## Virus Lineage : Bird Species
    real<lower=0> sigmasq_beta_Day_VLBS;
}

transformed parameters {

  ### Square root for convenience
    real<lower=0> sigma_alpha_IE;       
    real<lower=0> sigma_beta_Day_IE;
    real<lower=0> sigma_alpha_CIT; 
    real<lower=0> sigma_beta_Day_CIT;
    real<lower=0> sigma_alpha_VLBS;       
    real<lower=0> sigma_beta_Day_VLBS;
    real<lower=0> sigma_beta_Titer_VLBS;

    vector[N] lin_pred;   			 ## Linear predictor for survival
    
  ### 1 / Sqaure root
    sigma_alpha_IE = 1.0 / sqrt(sigmasq_alpha_IE);
    sigma_beta_Day_IE = 1.0 / sqrt(sigmasq_beta_Day_IE);
    sigma_alpha_CIT = 1.0 / sqrt(sigmasq_alpha_CIT);
    sigma_beta_Day_CIT = 1.0 / sqrt(sigmasq_beta_Day_CIT);
    sigma_alpha_VLBS = 1.0 / sqrt(sigmasq_alpha_VLBS);
    sigma_beta_Day_VLBS = 1.0 / sqrt(sigmasq_beta_Day_VLBS);
    
    for (i in 1:N) {
     lin_pred[i] =       
     alpha[VL[i]] +       		          ## intercept for each virus lineage
     alpha_IE_r[IE[i]] +             		  ## random intercept for each infection experiment
     alpha_CIT_r[CIT[i]] +           		  ## for each citation
     alpha_VLBS_r[VL_BS[i]] +        	          ## for each Bird Species
     (beta_Day[VL[i]] +           		  ## day slope for each virus 
     beta_Day_IE_r[IE[i]] +         		  ## random component of day slope, infection experiment
     beta_Day_VLBS_r[VL_BS[i]] +                  ## for Virus Lineage : Bird Species
     beta_Day_CIT_r[CIT[i]]) *                    ## for each citation
     Day[i] +            			  ## across days
     beta_LD[VL[i]] * LD[i];		          ## Log Dose
    }

}
model {

   alpha ~ normal(0.0, 1.0E3);		 	 ## diffuse priors for fixed effects
   beta_Day ~ normal(0.0, 1.0E3);
   beta_LD ~ normal(0.0, 1.0E3);

   sigmasq_alpha_IE ~ gamma(1.0E-3, 1.0E-3);	  ## random effects variances
   sigmasq_beta_Day_IE ~ gamma(1.0E-3, 1.0E-3);
   sigmasq_alpha_CIT ~ gamma(1.0E-3, 1.0E-3);
   sigmasq_beta_Day_CIT ~ gamma(1.0E-3, 1.0E-3);
   sigmasq_alpha_VLBS ~ gamma(1.0E-3, 1.0E-3);
   sigmasq_beta_Day_VLBS ~ gamma(1.0E-3, 1.0E-3);

  for (l in 1:N_IE) {
    alpha_IE_r[l] ~ normal(0, sigma_alpha_IE);
    beta_Day_IE_r[l] ~ normal(0, sigma_beta_Day_IE);
  }

  for (n in 1:N_CIT) {
    alpha_CIT_r[n] ~ normal(0, sigma_alpha_CIT);
    beta_Day_CIT_r[n] ~ normal(0, sigma_beta_Day_CIT);
  }

  for (q in 1:N_VL_BS) {
    alpha_VLBS_r[q] ~ normal(0, sigma_alpha_VLBS);
    beta_Day_VLBS_r[q] ~ normal(0, sigma_beta_Day_VLBS);
  }

  for (i in 1:N) {
     N_Ali[i] ~ binomial_logit(N_Samp[i], lin_pred[i]);
  }
   
}