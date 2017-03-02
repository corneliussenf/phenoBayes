data {
  int<lower=0> N; // Number of observations
  int<lower=0> NY; // Number of years
  int<lower=0> NP; // Number of pixels
  real doy[N]; // Day of year (predictor)
  real vi[N]; // Vegetation index (response)
  int<lower=0, upper=NY> year[N]; // Observation year
  int<lower=0, upper=NP> pixel[N]; // Observation pixel
  matrix[5, NP] mu_beta; # Centering of beta
  row_vector[5] sigma_beta_scale; # Scaling of sigma_beta
}

parameters {
  
  // Matrix for beta (scaled)
  matrix[5, NP] beta_raw;
  
  // Correlation matrix for beta with Cholesky factorization
  cholesky_factor_corr[5] cor_beta_L;
  
  // Scale parameters for beta (centered)
  vector<lower=0>[5] sigma_beta;
  
  // Prediction error scale
  real<lower=0> sigma;
  
  // Phi
  vector[NY] phi;
  
  // Scale parameters for phi with Cholesky factorization
  real<lower=0> sigma_phi;
  
}

transformed parameters {
  
  // Matrix for beta (non-centered)
  matrix[5, NP] beta;
  
  // Correlation matrix for beta
  matrix[5, 5] cor_beta;
  
  // Scaled sigma_beta
  vector[5] sigma_beta_scaled;
  
  // Scale sigma_beta
  for(i in 1:5)
    sigma_beta_scaled[i] = sigma_beta[1] * sigma_beta_scale[i];
    
  // Cholesky factorization of beta (incl. centering of beta)
  beta = diag_pre_multiply(sigma_beta_scaled, cor_beta_L) * beta_raw + mu_beta;
  
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
  sigma_beta ~ cauchy(0, 1);
  
  // Prior for phi
  phi ~ normal(0, sigma_phi);
  
  // Prior for scale of phi
  sigma_phi ~ normal(0, 5);
  
  // Prior for residual variance
  sigma ~ cauchy(0, 1);
  
  // Likelihood  
  for(i in 1:N)
    y[i] = (beta[1, pixel[i]]) +
        ((beta[2, pixel[i]] - (beta[5, pixel[i]] * doy[i])) / 
        (1 + exp(beta[3, pixel[i]] * (doy[i] - (beta[4, pixel[i]] + phi[year[i]])))));
  
  vi ~ normal(y, sigma);
  
}
