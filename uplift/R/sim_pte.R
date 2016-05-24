##################----------------------##################
# Simulations for personalized treatment effects
# "A Simple Method for Detecting Interactions between a
# Treatment and a Large Number of Covariates", Tian et al.
# (2012)
##################----------------------##################

# n = number of samples
# p = number of predictors
# ro = covariance between predictors 
# sigma = mutiplier of the error term
# beta.den = beta is mutiplied by 1/beta.den

sim_pte <- function(n = 1000, p = 20, rho = 0, sigma =  sqrt(2), beta.den = 4) {
  
  ### Check arguments
  if (n < 2)
    stop("uplift: The number of observations must be greater than 2")

  if (p < 4)
    stop("uplift: The number predictors must be equal or greater than 4")

  if (rho < 0 | rho > 1)
    stop("uplift: rho must be between 0 and 1")
  
  if (sigma < 0)
    stop("uplift: beta.den must be equal or greater than 0")
  
  if (beta.den <= 0)
    stop("uplift: beta.den must be greater than 0")
  
  ### Generate x ~ N~p(0, rho)
  cov.mat <- matrix(rho, nrow = p, ncol = p)
  cov.mat[row(cov.mat) == col(cov.mat)] <- 1
  x <- mvrnorm(n = n, mu = rep(0, p), Sigma = cov.mat)

  ### Random error from N~1(0,1)
  eps = rnorm(n, mean = 0, sd = 1)

  ### Treatment generated with equal probability at random 
  treat <- c(rep(-1, n/2), rep(1, n/2))
  treat <- treat[order(runif(n))]

  ### Main effects
  beta <- numeric(p)

  for (j in 1:p) {  
    beta[j] <- (-1)^j * (j >= 3 & j <= 10) / beta.den
  }

  ### Interaction effects
  gamma <- c(0.5, -0.5, 0.5, -0.5, rep(0, p - 4))

  ### Response variable
  y <- 1 * (x %*% beta + (x * treat) %*% gamma + sigma * eps > 0)

  ### "True" score
  ts <- pnorm(x %*%  (beta + gamma) / sigma) - pnorm(x %*%  (beta - gamma) / sigma)

  ### Returned value
  res <- data.frame(y, treat, x, ts = ts)

  return(res)

}

### END FUN


