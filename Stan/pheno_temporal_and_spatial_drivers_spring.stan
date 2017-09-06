data {
  int<lower=0> N; // Number of observations
  int<lower=0> NY; // Number of years
  int<lower=0> NP; // Number of pixels
  real doy[N]; // Day of year 
  real vi[N]; // Vegetation index 
  int<lower=0, upper=NY> year[N]; // Observation year
  int<lower=0, upper=NP> pixel[N]; // Observation pixel
  matrix[4, NP] mu_beta; // Centering of beta
  int lambda_null_mean; // Centering of lambda
  int<lower=0> PC_LS; // Number of landscape controls
  matrix[NP, PC_LS] landscape_controls; // Model matrix of landscape controls
  int<lower=0> PC_CLIM; // Number of climate controls
  matrix[NY, PC_CLIM] climate_controls; // Model matrix of climate controls (including intercept)
}

parameters {
  
  // Matrix for beta (scaled)
  matrix[4, NP] z_beta;
  
  // Correlation matrix for beta with Cholesky factorization
  cholesky_factor_corr[4] cor_beta_L;
  
  // Scale parameters for beta (centered)
  vector<lower=0>[4] sigma_beta;
  
  // Prediction error scale
  real<lower=0> sigma;
  
  // Phi (year-level hierachical model; centered)
  vector[NY] z_phi;
  
  // Scale parameters for phi
  real<lower=0> sigma_phi;
  
  // Climate control parameters
  vector[PC_CLIM] rho;
  
  // Start of season
  vector[NP] z_gamma;
  
  // Scale parameter for start of season
  real<lower=0> sigma_gamma;
  
  // Landscape control parameters
  vector[PC_LS] lambda;
  real lambda_null;
  
}

transformed parameters {
  
  // Matrix for beta
  matrix[4, NP] beta;
  
  // Correlation matrix for beta
  matrix[4, 4] cor_beta;
  
  // Start of season
  vector[NP] gamma;
  
  // Gamma mean estimate (landscape controls)
  vector[NP] landscape_model;
  
  // Vector for phi (non-centered)
  vector[NY] phi;
  
  // Phi mean estimate (climate controls)
  vector[NY] climate_model;
  
  // Cholesky factorization of beta (incl. centering of beta)
  beta = diag_pre_multiply(sigma_beta * 0.1, cor_beta_L) * z_beta + mu_beta;
  
  // Backtransform Cholesky correlation matrix (LL') of beta 
  cor_beta = tcrossprod(cor_beta_L);
  
  // Model for landscape controls
  landscape_model = landscape_controls * lambda + (lambda_null + lambda_null_mean);
  
  // Reparameterize gamma
  gamma = landscape_model + z_gamma * sigma_gamma;
  
  // Model for phi (climate controls)
  climate_model = climate_controls * rho;
  
  // Reparameterize phi
  phi = climate_model + z_phi * sigma_phi;
  
}

model {
  
  vector[N] y;
  
  // Prior for beta
  to_vector(z_beta) ~ normal(0, 1);
  
  // Prior for correlation matrix of beta
  cor_beta_L ~ lkj_corr_cholesky(2.0);
  
  // Prior for scales of beta
  sigma_beta ~ cauchy(0, 1);
  
  // Priors for landscape controls
  lambda ~ student_t(3, 0, 1);
  
  // Prior on start of season (scaled)
  lambda_null ~ normal(0, 1);
  
  // Priors for climate controls
  rho ~ student_t(3, 0, 1);
  
  // Prior for gamma
  z_gamma ~ normal(0, 1);
  
  // Prior for scale of thetha
  sigma_gamma ~ cauchy(0, 5);
  
  // Prior for phi
  z_phi ~ normal(0, 1);
  
  // Prior for scale of phi
  sigma_phi ~ cauchy(0, 5);
  
  // Prior for residual variance
  sigma ~ cauchy(0, 1);
  
  // Likelihood  
  for(i in 1:N)
    y[i] = (beta[1, pixel[i]]) + ((beta[2, pixel[i]]) - (beta[4, pixel[i]]) * doy[i]) * 
        (1 / (1 + exp(-(beta[3, pixel[i]]) * (doy[i] - (gamma[pixel[i]] + phi[year[i]])))));
  
  vi ~ normal(y, sigma);
  
}

generated quantities {
  
  vector[N] sim;
  
  # Posterior simulations
  for (n in 1:N) 
    sim[n] = normal_rng((beta[1, pixel[n]]) + ((beta[2, pixel[n]]) - (beta[4, pixel[n]]) * doy[n]) * 
          (1 / (1 + exp(-(beta[3, pixel[n]]) * (doy[n] - (gamma[pixel[n]] + phi[year[n]]))))), sigma);
  
}

