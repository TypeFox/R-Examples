
##################################
## models for latent parameters ##
##################################

# Generate AR(1) variates

r_AR1 <- function(iterations, series_length, rho, sigma_sq) 
  matrix(rnorm(iterations * series_length), iterations, series_length) %*% chol(sigma_sq * rho^as.matrix(dist(1:series_length)))


# solve for conditional prevalence 

logit <- function(x) log(x) - log(1 - x)
expit <- function(x) 1 / (1 + exp(-x))

E_logitnorm <- function(eta_star, sigma_sq)
  sapply(eta_star, function(eta_star)
    integrate(function(v) 
      expit(eta_star + v) * exp(-v^2 / (2 * sigma_sq)) / sqrt(2 * pi * sigma_sq), 
              lower = -Inf, upper = Inf)$value)

eta_star_phi <- function(marginal_mean, sigma_sq, interval = c(-1000,1000)) {
  n <- max(length(marginal_mean), length(sigma_sq))
  marginal_mean <- rep(marginal_mean, length.out = n)
  sigma_sq <- rep(sigma_sq, length.out = n)
  mapply(function(marginal_mean, sigma_sq)
    uniroot(function(x) E_logitnorm(x, sigma_sq) - marginal_mean, 
            interval = interval)$root,
         marginal_mean = marginal_mean, sigma_sq = sigma_sq)
}


# generate random, dependent phi values

r_phi_star <- function(iterations, series_length, phi_marg, rho, sigma_sq) {
  eta_cond <- rep(eta_star_phi(phi_marg, sigma_sq), length.out = series_length)
  nu <- t(r_AR1(iterations, series_length, rho, sigma_sq))
  expit(eta_cond + nu)
}

# generate random, dependent zeta values

r_zeta_star <- function(iterations, series_length, zeta_marg, rho, sigma_sq) {
  nu <- t(r_AR1(iterations, series_length, rho, sigma_sq)) - sigma_sq / 2
  zeta_marg * exp(nu)
}

# smooth covariance matrix
smooth_cov <- function(V) {
  n <- dim(V)[1]
  smooth_cov <- sapply(1:n, function(x) ifelse(x < n, mean(diag(V[x:n,1:(n-x+1)])), V[x,1]))
  matrix(smooth_cov[as.matrix(dist(1:n, diag=TRUE)) + 1], n, n) / smooth_cov[1]
}

