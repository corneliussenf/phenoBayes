data {
  int<lower=0> N; // Number of observations
  int<lower=0> NY; // Number of years
  int<lower=0> NP; // Number of pixels
  real doy[N]; // Day of year
  real vi[N]; // Vegetation index
  int<lower=0, upper=NY> year[N]; // Observation year (consecutive number 1, ..., NY)
  int<lower=0, upper=NP> pixel[N]; // Observation pixel (consecutive number 1, ..., NP)
  matrix[5, NP] mu_beta; # Centering of beta
  row_vector[5] sigma_beta_scale; # Scaling of sigma_beta
}

parameters {
  
  // Matrix for beta (centered)
  matrix[5, NP] beta_raw;
  
  // Correlation matrix for beta with Cholesky factorization
  cholesky_factor_corr[5] cor_beta_L;
  
  // Scale parameters for beta (non-scaled)
  vector<lower=0>[5] sigma_beta_raw;
  
  // Prediction error scale
  real<lower=0> sigma;
  
  // Vector for phi
  vector[NY] phi;
  
  // Scale parameters for phi
  real<lower=0> sigma_phi;
  
}

transformed parameters {
  
  // Matrix for beta (non-centered)
  matrix[5, NP] beta;
  
  // Correlation matrix for beta
  matrix[5, 5] cor_beta;
  
  // Sigma_beta (scaled)
  vector[5] sigma_beta;
  
  // Scale sigma_beta
  for(i in 1:5)
    sigma_beta[i] = sigma_beta_raw[i] * sigma_beta_scale[i];
    
  // Cholesky factorization of beta (incl. centering of beta)
  beta = diag_pre_multiply(sigma_beta, cor_beta_L) * beta_raw + mu_beta;
  
  // Backtransform Cholesky correlation matrix (LL') of beta 
  cor_beta = tcrossprod(cor_beta_L);
  
}

model {
  
  vector[N] y;
  
  // Prior for beta
  to_vector(beta_raw) ~ normal(0, 1);
  
  // Prior for correlation matrix of beta
  cor_beta_L ~ lkj_corr_cholesky(2.0);
  
  // Prior for scales of beta
  sigma_beta_raw ~ cauchy(0, 1);
  
  // Prior for phi
  phi ~ normal(0, sigma_phi);
  
  // Prior for scale of phi
  sigma_phi ~ normal(0, 5);
  
  // Prior prediction error scale
  sigma ~ cauchy(0, 1);
  
  // Likelihood  
  for(i in 1:N)
    y[i] = (beta[1, pixel[i]]) + ((beta[2, pixel[i]]) - (beta[5, pixel[i]]) * doy[i]) * 
        (1 / (1 + exp(-(beta[3, pixel[i]]) * (doy[i] - (beta[4, pixel[i]] + phi[year[i]])))));
  
  vi ~ normal(y, sigma);
  
}
