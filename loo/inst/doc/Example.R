## ---- eval=FALSE---------------------------------------------------------
#  'logistic.stan'
#  
#  data {
#    int<lower=0> N; # number of data points
#    int<lower=0> P; # number of predictors (including intercept)
#    int<lower=0,upper=1> y[N]; # binary outcome
#    matrix[N,P] x; # predictors (including intercept)
#    real a;
#  }
#  parameters {
#    real beta0;
#    vector[P] beta;
#  }
#  model {
#    beta0 ~ student_t(7, a, 0.1);
#    beta ~ student_t(7, 0, 1);
#    y ~ bernoulli_logit(beta0 + x * beta);
#  }
#  generated quantities {
#    vector[N] log_lik;
#    for (n in 1:N)
#      log_lik[n] <- bernoulli_logit_log(y[n], beta0 + x[n] * beta);
#  }

## ---- eval=FALSE---------------------------------------------------------
#  library("rstan")
#  
#  # Prepare data
#  url <- "http://stat.columbia.edu/~gelman/arm/examples/arsenic/wells.dat"
#  wells <- read.table(url)
#  wells$dist100 <- wells$dist / 100 # rescale
#  y <- wells$switch
#  a <- qlogis(mean(y)) # i.e., a = logit(Pr(y = 1))
#  x <- scale(model.matrix(~ 0 + dist + arsenic, wells))
#  data <- list(N = nrow(x), P = ncol(x), a = a, x = x, y = y)
#  
#  # Fit model
#  fit1 <- stan("logistic.stan", data = data)

## ---- eval=FALSE---------------------------------------------------------
#  library("loo")
#  
#  # Extract log-likelihood and compute LOO
#  log_lik1 <- extract_log_lik(fit1)
#  loo1 <- loo(log_lik1) # or waic(log_lik1) to compute WAIC

## ---- eval=FALSE---------------------------------------------------------
#  # First run a second model using log(arsenic) instead of arsenic
#  data$x <- scale(model.matrix(~ 0 + dist100 + log(arsenic), wells))
#  fit2 <- stan(fit = fit1, data = data)
#  log_lik2 <- extract_log_lik(fit2)
#  loo2 <- loo(log_lik2)
#  
#  # Compare
#  diff <- compare(loo1, loo2)

