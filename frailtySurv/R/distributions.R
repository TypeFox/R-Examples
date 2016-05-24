# This file contains R implementations of frailty distribution functions
# The corresponding C++ functions are much faster and should be used instead
# Any function ending with _r has a corresponding c function ending with _c
# The function here are provided for mainly for testing purposes. Most functions
# have a _numeric counterpart, where results are obtained numerically.

################################################################################

# Example of calculating phi_ with numerical integrals
phi_numerical <- function(k, N_dot, H_dot, density_params, density_fn) {
  f <- function(w) {
    w^(N_dot + k - 1) * exp(-w*H_dot) * density_fn(w, density_params)
  }
  integrate(f, 0, Inf)
}

# Example of calculating phi_ with the LT
phi_laplace <- function(k, N_dot, H_dot, density_params, density_LT) {
  density_LT(N_dot + k - 1, H_dot, density_params)*(-1)^(N_dot + k - 1)
}

# Laplace transform, evaluated numerically
lt_numeric <- function(m, s, theta, density_fn) {
  integrate(function(t) {
    (-t)^m * exp(-s * t) * density_fn(t, theta)
  }, 0, Inf)$value
}

################################################################################
# Gamma distribution

# 
# Gamma density
# 
dgamma_r <- function(x, theta) {
  dgamma(x, 1/theta, 1/theta)
}

# 
# Gamma random variable generation
# By LT: rlaptrans(n, lt_dgamma_r, p=0, theta=theta)
# 
rgamma_r <- function(n, theta) {
  if (theta == 0)
    rep(1, n)
  else
    rgamma(n, 1/theta, 1/theta)
}

#
# Gamma variance
#
vgamma_r <- function(theta) {
  theta
}

# 
# Deriv of gamma density wrt. theta
# 
deriv_dgamma_r <- function(x, theta) {
  # deriv_idx ignored here since there is only one parameter
  ( (x/theta)^(1/theta - 1) * 
      exp(-x/theta) * 
      (log(theta/x) + digamma(1/theta) + x - 1) )/
    (gamma(1/theta)*theta^3)
}

# 
# Deriv of gamma density wrt. theta, evaluated numerically
# 
deriv_dgamma_r_numeric <- function(x, theta) {
  grad(function(theta) dgamma_r(x, theta), theta)
}

# 
# 2nd Deriv of gamma density wrt. theta, evaluated numerically
# 
deriv_deriv_dgamma_r <- function(x, theta) {
  (((x/theta)^(1/theta - 1) * (exp(-x/theta) * (x/theta^2)) - ((x/theta)^(1/theta - 
  1) * (log((x/theta)) * (1/theta^2)) + (x/theta)^((1/theta - 
  1) - 1) * ((1/theta - 1) * (x/theta^2))) * exp(-x/theta)) * 
   (log(theta/x) + digamma(1/theta) + x - 1) + (x/theta)^(1/theta - 
  1) * exp(-x/theta) * (1/x/(theta/x) - 1/theta^2 * trigamma(1/theta)))/(gamma(1/theta) * 
   theta^3) - ((x/theta)^(1/theta - 1) * exp(-x/theta) * (log(theta/x) + 
  digamma(1/theta) + x - 1)) * (gamma(1/theta) * (3 * theta^2) - 
  1/theta^2 * (gamma(1/theta) * digamma(1/theta)) * theta^3)/(gamma(1/theta) * 
  theta^3)^2
}

#
# Laplace transform of the gamma density, as defined above
# 
lt_dgamma_r <- function(m, s, theta) {
  (-1)^m * 
    (1/theta)^(1/theta) * 
    ((1/theta) + s)^(-((1/theta) + m)) * 
    gamma((1/theta) + m)/gamma((1/theta))
}

#
# Derivative of the Laplace transform of the gamma density, wrt theta,
# pth derivative wrt. s
# 
deriv_lt_dgamma_r <- function(m, s, theta) {
  (((-1)^m * (1/theta)^(1/theta) * (((1/theta) + s)^(-((1/theta) + 
  m)) * (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + 
  s)^((-((1/theta) + m)) - 1) * ((-((1/theta) + m)) * (1/theta^2))) - 
  (-1)^m * ((1/theta)^(1/theta) * (log((1/theta)) * (1/theta^2)) + 
  (1/theta)^((1/theta) - 1) * ((1/theta) * (1/theta^2))) * 
  ((1/theta) + s)^(-((1/theta) + m))) * gamma((1/theta) + 
  m) - (-1)^m * (1/theta)^(1/theta) * ((1/theta) + s)^(-((1/theta) + 
  m)) * (1/theta^2 * (gamma((1/theta) + m) * digamma((1/theta) + 
  m))))/gamma((1/theta)) + (-1)^m * (1/theta)^(1/theta) * ((1/theta) + 
  s)^(-((1/theta) + m)) * gamma((1/theta) + m) * (1/theta^2 * 
  (gamma((1/theta)) * digamma((1/theta))))/gamma((1/theta))^2
}

# 
# Gamma LT 2nd derivative wrt. theta. This obviously was not solved by hand.
# 
deriv_deriv_lt_dgamma_r <- function(m, s, theta) {
  (((-1)^m * (1/theta)^(1/theta) * ((((1/theta) + s)^(-((1/theta) + 
  m)) * (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + 
  s)^((-((1/theta) + m)) - 1) * ((-((1/theta) + m)) * (1/theta^2))) * 
  (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + s)^(-((1/theta) + 
  m)) * (log(((1/theta) + s)) * (2 * theta/(theta^2)^2) + 1/theta^2/((1/theta) + 
   s) * (1/theta^2)) - ((((1/theta) + s)^((-((1/theta) + m)) - 
  1) * (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + 
   s)^(((-((1/theta) + m)) - 1) - 1) * (((-((1/theta) + m)) - 
   1) * (1/theta^2))) * ((-((1/theta) + m)) * (1/theta^2)) + 
  ((1/theta) + s)^((-((1/theta) + m)) - 1) * (1/theta^2 * (1/theta^2) - 
  (-((1/theta) + m)) * (2 * theta/(theta^2)^2)))) - (-1)^m * 
  ((1/theta)^(1/theta) * (log((1/theta)) * (1/theta^2)) + (1/theta)^((1/theta) - 
   1) * ((1/theta) * (1/theta^2))) * (((1/theta) + s)^(-((1/theta) + 
   m)) * (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + 
   s)^((-((1/theta) + m)) - 1) * ((-((1/theta) + m)) * (1/theta^2))) - 
  ((-1)^m * ((1/theta)^(1/theta) * (log((1/theta)) * (1/theta^2)) + 
   (1/theta)^((1/theta) - 1) * ((1/theta) * (1/theta^2))) * 
   (((1/theta) + s)^(-((1/theta) + m)) * (log(((1/theta) + 
   s)) * (1/theta^2)) - ((1/theta) + s)^((-((1/theta) + 
  m)) - 1) * ((-((1/theta) + m)) * (1/theta^2))) - 
   (-1)^m * ((1/theta)^((1/theta) - 1) * ((1/theta) * (2 * 
   theta/(theta^2)^2) + 1/theta^2 * (1/theta^2)) + ((1/theta)^((1/theta) - 
   1) * (log((1/theta)) * (1/theta^2)) + (1/theta)^(((1/theta) - 
   1) - 1) * (((1/theta) - 1) * (1/theta^2))) * ((1/theta) * 
   (1/theta^2)) + ((1/theta)^(1/theta) * (log((1/theta)) * 
  (2 * theta/(theta^2)^2) + 1/theta^2/(1/theta) * (1/theta^2)) + 
   ((1/theta)^(1/theta) * (log((1/theta)) * (1/theta^2)) + 
  (1/theta)^((1/theta) - 1) * ((1/theta) * (1/theta^2))) * 
   (log((1/theta)) * (1/theta^2)))) * ((1/theta) + 
   s)^(-((1/theta) + m)))) * gamma((1/theta) + m) - 
   ((-1)^m * (1/theta)^(1/theta) * (((1/theta) + s)^(-((1/theta) + 
   m)) * (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + 
   s)^((-((1/theta) + m)) - 1) * ((-((1/theta) + m)) * (1/theta^2))) - 
  (-1)^m * ((1/theta)^(1/theta) * (log((1/theta)) * (1/theta^2)) + 
  (1/theta)^((1/theta) - 1) * ((1/theta) * (1/theta^2))) * 
  ((1/theta) + s)^(-((1/theta) + m))) * (1/theta^2 * 
   (gamma((1/theta) + m) * digamma((1/theta) + m))) - (((-1)^m * 
  (1/theta)^(1/theta) * (((1/theta) + s)^(-((1/theta) + m)) * 
   (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + s)^((-((1/theta) + 
  m)) - 1) * ((-((1/theta) + m)) * (1/theta^2))) - (-1)^m * 
  ((1/theta)^(1/theta) * (log((1/theta)) * (1/theta^2)) + (1/theta)^((1/theta) - 
   1) * ((1/theta) * (1/theta^2))) * ((1/theta) + s)^(-((1/theta) + 
  m))) * (1/theta^2 * (gamma((1/theta) + m) * digamma((1/theta) + 
  m))) - (-1)^m * (1/theta)^(1/theta) * ((1/theta) + s)^(-((1/theta) + 
   m)) * (1/theta^2 * (gamma((1/theta) + m) * (1/theta^2 * trigamma((1/theta) + 
  m)) + 1/theta^2 * (gamma((1/theta) + m) * digamma((1/theta) + 
  m)) * digamma((1/theta) + m)) + 2 * theta/(theta^2)^2 * (gamma((1/theta) + 
   m) * digamma((1/theta) + m)))))/gamma((1/theta)) + (((-1)^m * 
  (1/theta)^(1/theta) * (((1/theta) + s)^(-((1/theta) + m)) * 
   (log(((1/theta) + s)) * (1/theta^2)) - ((1/theta) + s)^((-((1/theta) + 
  m)) - 1) * ((-((1/theta) + m)) * (1/theta^2))) - (-1)^m * 
  ((1/theta)^(1/theta) * (log((1/theta)) * (1/theta^2)) + (1/theta)^((1/theta) - 
   1) * ((1/theta) * (1/theta^2))) * ((1/theta) + s)^(-((1/theta) + 
  m))) * gamma((1/theta) + m) - (-1)^m * (1/theta)^(1/theta) * 
   ((1/theta) + s)^(-((1/theta) + m)) * (1/theta^2 * (gamma((1/theta) + 
  m) * digamma((1/theta) + m)))) * (1/theta^2 * (gamma((1/theta)) * 
   digamma((1/theta))))/gamma((1/theta))^2 + (((((-1)^m * (1/theta)^(1/theta) * 
   (((1/theta) + s)^(-((1/theta) + m)) * (log(((1/theta) + s)) * 
  (1/theta^2)) - ((1/theta) + s)^((-((1/theta) + m)) - 
  1) * ((-((1/theta) + m)) * (1/theta^2))) - (-1)^m * ((1/theta)^(1/theta) * 
   (log((1/theta)) * (1/theta^2)) + (1/theta)^((1/theta) - 1) * 
   ((1/theta) * (1/theta^2))) * ((1/theta) + s)^(-((1/theta) + 
   m))) * gamma((1/theta) + m) - (-1)^m * (1/theta)^(1/theta) * 
  ((1/theta) + s)^(-((1/theta) + m)) * (1/theta^2 * (gamma((1/theta) + 
   m) * digamma((1/theta) + m)))) * (1/theta^2 * (gamma((1/theta)) * 
  digamma((1/theta)))) - (-1)^m * (1/theta)^(1/theta) * ((1/theta) + 
   s)^(-((1/theta) + m)) * gamma((1/theta) + m) * (1/theta^2 * 
   (gamma((1/theta)) * (1/theta^2 * trigamma((1/theta))) + 1/theta^2 * 
  (gamma((1/theta)) * digamma((1/theta))) * digamma((1/theta))) + 
   2 * theta/(theta^2)^2 * (gamma((1/theta)) * digamma((1/theta)))))/gamma((1/theta))^2 + 
  (-1)^m * (1/theta)^(1/theta) * ((1/theta) + s)^(-((1/theta) + 
  m)) * gamma((1/theta) + m) * (1/theta^2 * (gamma((1/theta)) * 
   digamma((1/theta)))) * (2 * (1/theta^2 * (gamma((1/theta)) * 
   digamma((1/theta))) * gamma((1/theta))))/(gamma((1/theta))^2)^2)
}

# 
# Gamma LT derivative wrt. theta, evaluated numerically
# 
deriv_lt_dgamma_r_numeric <- function(m, s, theta) {
  grad(function(theta) lt_dgamma_r(m, s, theta), theta)
}

# 
# Gamma LT 2nd derivative wrt. theta, evaluated numerically
# 
deriv_deriv_lt_dgamma_r_numeric <- function(m, s, theta) {
  grad(function(theta) deriv_lt_dgamma_r_numeric(m, s, theta), theta)
}

# 
# pth moment of the gamma distribution, evaluated by the gamma density LT
# 
lt_dgamma_moment_r <- function(p, theta) {
  (-1)^p * lt_dgamma_r(p, 0, theta)
}

################################################################################
# Positive stable distribution. Density and random generation only

# 
# Positive stable density]
# 
dposstab_r <- function(x, alpha, K=100) {
  
  f <- function(x) {
    k <- 1:K
    -1 / (pi * x) * sum(
      gamma(k * alpha + 1) / factorial(k) * 
        (-x ^ (-alpha)) ^ k * sin(alpha * k * pi))
  }
  
  Vectorize(f)(x)
}

# 
# Random generation from a positive stable distribution
#
rposstab_r <- function(n, alpha) {
  if (alpha == 1)
    rep(1, n)
  else
    rlaptrans(n, lt_dposstab_r, p=0, alpha=alpha)
}

#
# Laplace transform of the positive stable distribution
# 
lt_dposstab_r <- function(p, s, alpha) {
  if (p == 0) {
    return(exp(-alpha*(s^alpha)/alpha))
  }
  
  (-1)^p * lt_dposstab_r(0, s, alpha) * sum(vapply(1:p, function(j)
    lt_dpvf_coef_r(p, j, alpha) * alpha^j * s^(j*alpha - p), 0
  ))
}

################################################################################
# Power variance function distribution (PVF)

#
# PVF density
# 
dpvf_r <- function(x, alpha, K=100) {
  Vectorize(function(x) {
    dposstab_r(x*alpha^(1/alpha), alpha, K)*alpha^(1/alpha)*exp(-x)*exp(1/alpha)
  })(x)
}

# 
# Generate samples from a PVF using its LT
# 
rpvf_r <- function(n, alpha) {
  if (alpha == 1)
    rep(1, n)
  else
    rlaptrans(n, lt_dpvf_r, p=0, alpha=alpha)
}

# 
# PVF variance
# 
vpvf_r <- function(alpha) {
  1 - alpha
}

#
# Coefficients for the PVF LT derivatives
# 
lt_dpvf_coef_r <- function(p, j, alpha) {
  if (p == j) return(1)
  if (j == 1) return(gamma(p - alpha)/gamma(1 - alpha))
  
  return(lt_dpvf_coef_r(p - 1, j - 1, alpha) + 
           lt_dpvf_coef_r(p - 1, j, alpha)*((p - 1) - j*alpha))
}

# 
# Coefficients for the PVF LT
# 
deriv_lt_dpvf_coef_r <- function(p, j, alpha) {
  
  if (p == j) return(0)
  
  if (j == 1) {
    numer <- gamma(p - alpha)*(psigamma(1 - alpha) - psigamma(p - alpha))
    denom <- gamma(1 - alpha)
    return(numer/denom)
  }
  
  return(deriv_lt_dpvf_coef_r(p - 1, j - 1, alpha) + 
           deriv_lt_dpvf_coef_r(p - 1, j, alpha) * ((p - 1) - j*alpha) -
           j * lt_dpvf_coef_r(p - 1, j, alpha)) # last coef is not a deriv coef
}

# 
# Coefficients for the PVF LT, numerical deriv wrt. alpha
# 
deriv_lt_dpvf_coef_r_numeric <- function(p, j, alpha) {
  return(grad(function(alpha) lt_dpvf_coef(p, j, alpha), alpha))
}

# 
# Coefficients for the PVF LT
# 
deriv_deriv_lt_dpvf_coef_r <- function(m, j, alpha) {

  if (m == j) return(0)
  
  if (j == 1) {
    return((gamma(m - alpha) * trigamma(m - alpha) + gamma(m - alpha) * 
    digamma(m - alpha) * digamma(m - alpha))/gamma(1 - alpha) - 
     gamma(m - alpha) * digamma(m - alpha) * (gamma(1 - alpha) * 
    digamma(1 - alpha))/gamma(1 - alpha)^2 - ((gamma(m - 
     alpha) * (gamma(1 - alpha) * trigamma(1 - alpha) + gamma(1 - 
    alpha) * digamma(1 - alpha) * digamma(1 - alpha)) + gamma(m - 
    alpha) * digamma(m - alpha) * (gamma(1 - alpha) * digamma(1 - 
    alpha)))/gamma(1 - alpha)^2 - gamma(m - alpha) * (gamma(1 - 
    alpha) * digamma(1 - alpha)) * (2 * (gamma(1 - alpha) * digamma(1 - 
    alpha) * gamma(1 - alpha)))/(gamma(1 - alpha)^2)^2))
}
  
  term1 <- deriv_deriv_lt_dpvf_coef_r(m - 1, j - 1, alpha)
  term2 <- deriv_deriv_lt_dpvf_coef_r(m - 1, j, alpha) * ((m - 1) - j*alpha)
  term3 <- -j * deriv_lt_dpvf_coef_r(m - 1, j, alpha)
  
  term1 + term2 + 2*term3
}

# 
# Coefficients for the PVF LT, numerical deriv wrt. alpha
# 
deriv_deriv_lt_dpvf_coef_r_numeric <- function(p, j, alpha) {
  return(grad(function(alpha) deriv_lt_dpvf_coef_r(p, j, alpha), alpha))
}

# 
# Laplace transform of the one-parameter PVF distribution defined above
# 
lt_dpvf_r <- function(p, s, alpha) {
  if (p == 0) {
    return(exp(-((1 + s)^alpha - 1)/alpha))
  }
  
  (-1)^p * lt_dpvf_r(0, s, alpha) * sum(vapply(1:p, function(j)
    lt_dpvf_coef_r(p, j, alpha) * (1 + s)^(j*alpha - p), 0
  ))
}

# 
# pth moment of the one-parameter PVF distribution, as defined above
# 
lt_dpvf_moment_r <- function(alpha, p) { 
  (-1)^p * lt_dpvf_r(p, 0, alpha)
}

#
# Derivative wrt. alpha of the PVF LT, evaluated numerically
# 
deriv_lt_dpvf_r_numeric <- function(p, s, alpha) {
  grad(function(alpha) lt_dpvf_r(p, s, alpha), alpha)
}

# 
# PVF LT deriv wrt. alpha
# 
deriv_lt_dpvf_r <- function(m, s, alpha) {
  if (m == 0) {
    term1 <- ((s + 1)^alpha - 1)/alpha^2
    term2 <- ((s + 1)^alpha * log(1 + s))/alpha
    return(lt_dpvf_r(0, s, alpha) * (term1 - term2))
  }
  
  term1 <- (-1)^m * deriv_lt_dpvf_r(0, s, alpha) * 
    sum(vapply(1:m, function(j)
      lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m), 0))
  
  term2 <- (-1)^m * lt_dpvf_r(0, s, alpha) *
    sum(vapply(1:m, function(j)
      deriv_lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m) + 
        lt_dpvf_coef_r(m, j, alpha) * j * (1 + s)^(j*alpha - m) * log(1 + s), 0))
  
  term1 + term2
}

# 
# PVF LT 2nd deriv wrt. alpha
# 
deriv_deriv_lt_dpvf_r <- function(m, s, alpha) {
  if (m == 0) {
    return(-(exp(-((1 + s)^alpha - 1)/alpha) * ((1 + s)^alpha * log((1 + 
     s)) * log((1 + s))/alpha - (1 + s)^alpha * log((1 + s))/alpha^2 - 
    ((1 + s)^alpha * log((1 + s))/alpha^2 - ((1 + s)^alpha - 
     1) * (2 * alpha)/(alpha^2)^2)) - exp(-((1 + s)^alpha - 
    1)/alpha) * ((1 + s)^alpha * log((1 + s))/alpha - ((1 + s)^alpha - 
     1)/alpha^2) * ((1 + s)^alpha * log((1 + s))/alpha - ((1 + 
     s)^alpha - 1)/alpha^2)))
  }
  
  lt <- lt_dpvf_r(0, s, alpha)
  dlt <- deriv_lt_dpvf_r(0, s, alpha)
  ddlt <- deriv_deriv_lt_dpvf_r(0, s, alpha)
  
  sum1 <- sum(vapply(1:m, function(j)
    lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m) * j^2 * log(1 + s)^2
    , 0))
  
  sum2 <- sum(vapply(1:m, function(j)
    deriv_lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m) * j * log(1 + s)
    , 0))
  
  sum3 <- sum(vapply(1:m, function(j)
    deriv_deriv_lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m)
    , 0))
  
  sum4 <- sum(vapply(1:m, function(j)
    lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m) * j * log(1 + s)
    , 0))
  
  sum5 <- sum(vapply(1:m, function(j)
    lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m)
    , 0))
  
  sum6 <- sum(vapply(1:m, function(j)
    deriv_lt_dpvf_coef_r(m, j, alpha) * (1 + s)^(j*alpha - m)
    , 0))
  
  term1 <- (-1)^m * lt * sum1
  term2 <- 2 * (-1)^m * lt * sum2
  term3 <- (-1)^m * lt * sum3
  term4 <- 2 * (-1)^m * dlt * sum4
  term5 <- (-1)^m * ddlt * sum5
  term6 <- 2 * (-1)^m * dlt * sum6
  
  term1 + term2 + term3 + term4 + term5 + term6
}

#
# PVF LT 2nd deriv wrt. alpha, evaluated numerically
# 
deriv_deriv_lt_dpvf_r_numeric <- function(p, s, alpha) {
  grad(function(alpha) deriv_lt_dpvf_r(p, s, alpha), alpha)
}

################################################################################
# Log-normal

# 
# Log-normal density
# 
dlognormal_r <- function(x, theta) {
  dlnorm(x, 0, sqrt(theta))
}

# 
# Log-normal random generation
# 
rlognormal_r <- function(n, theta) {
  if (theta == 0)
    rep(1, n)
  else
    rlnorm(n, 0, sqrt(theta))
}

# 
# Log-normal variance
# 
vlognormal_r <- function(theta) {
  exp(2*theta) - exp(theta)
}

# 
# Deriv Log-normal
# 
deriv_dlognormal_r <- function(x, theta) {
  term1_numer <- log(x)^2 * exp(-(log(x)^2)/(2*theta))
  term1_denom <- 2 * sqrt(2*pi) * theta^(5/2) * x
  term2_numer <- exp(-(log(x)^2)/(2*theta))
  term2_denom <- 2 * sqrt(2*pi) * theta^(3/2) * x
  term1_numer/term1_denom - term2_numer/term2_denom
}

# 
# Deriv log-normal, numerical
# 
deriv_dlognormal_r_numeric <- function(x, theta) {
  grad(function(theta) dlognormal_r(x, theta), theta)
}

# 
# 2nd deriv log-normal
# 
deriv_deriv_dlognormal_r <- function(x, theta) {
  log(x)^2 * (exp(-(log(x)^2)/(2 * theta)) * ((log(x)^2) * 2/(2 * 
  theta)^2))/(2 * sqrt(2 * pi) * theta^(5/2) * x) - (log(x)^2 * 
   exp(-(log(x)^2)/(2 * theta))) * (2 * sqrt(2 * pi) * (theta^((5/2) - 
   1) * (5/2)) * x)/(2 * sqrt(2 * pi) * theta^(5/2) * x)^2 - 
  (exp(-(log(x)^2)/(2 * theta)) * ((log(x)^2) * 2/(2 * theta)^2)/(2 * 
  sqrt(2 * pi) * theta^(3/2) * x) - (exp(-(log(x)^2)/(2 * 
  theta))) * (2 * sqrt(2 * pi) * (theta^((3/2) - 1) * (3/2)) * 
  x)/(2 * sqrt(2 * pi) * theta^(3/2) * x)^2)
}

# 
# 2nd deriv log-normal, numerical
# 
deriv_deriv_dlognormal_r_numeric <- function(x, theta) {
  grad(function(theta) deriv_dlognormal_r(x, theta), theta)
}

################################################################################
# Inverse Gaussian

dinvgauss_r <- function(x, theta) {
  exp(-((x - 1)^2)/(2*theta*x))/sqrt(2*pi*theta*(x^3))
}

rinvgauss_r <- function(n, theta) {
  if (theta == 0)
    rep(1, n)
  else {
    mean <- 1
    shape <- 1/theta
    dispersion <- 1 
    # From statmod::rinvgauss
    if (!is.null(shape)) 
      dispersion <- 1/shape
    if (length(n) > 1L) 
      n <- length(n)
    else n <- as.integer(n)
    if (n < 0L) 
      stop("n can't be negative")
    if (n == 0L || length(mean) == 0L || length(dispersion) == 
        0L) 
      return(numeric(0L))
    mu <- rep_len(mean, n)
    phi <- rep_len(dispersion, n)
    r <- rep_len(0, n)
    i <- (mu > 0 & phi > 0)
    i[is.na(i)] <- FALSE
    if (!all(i)) {
      r[!i] <- NA
      n <- sum(i)
    }
    phi[i] <- phi[i] * mu[i]
    Y <- rchisq(n, df = 1)
    X1 <- 1 + phi[i]/2 * (Y - sqrt(4 * Y/phi[i] + Y^2))
    firstroot <- as.logical(rbinom(n, size = 1L, prob = 1/(1 + 
                                                             X1)))
    r[i][firstroot] <- X1[firstroot]
    r[i][!firstroot] <- 1/X1[!firstroot]
    mu * r
  }
}

vinvgauss_r <- function(theta) {
  theta
}

deriv_dinvgauss_r_numeric <- function(x, theta) {
  grad(function(theta) dinvgauss_r(x, theta), theta)
}

deriv_dinvgauss_r <- function(x, theta) {
  exp(-((x - 1)^2)/(2 * theta * x)) * (((x - 1)^2) * (2 * x)/(2 * 
  theta * x)^2)/sqrt(2 * pi * theta * (x^3)) - exp(-((x - 1)^2)/(2 * 
   theta * x)) * (0.5 * (2 * pi * (x^3) * (2 * pi * theta * 
   (x^3))^-0.5))/sqrt(2 * pi * theta * (x^3))^2
}

deriv_deriv_dinvgauss_r_numeric <- function(x, theta) {
  grad(function(theta) deriv_dinvgauss_r(x, theta), theta)
}

deriv_deriv_dinvgauss_r <- function(x, theta) {
  (exp(-((x - 1)^2)/(2 * theta * x)) * (((x - 1)^2) * (2 * x)/(2 * 
   theta * x)^2) * (((x - 1)^2) * (2 * x)/(2 * theta * x)^2) - 
   exp(-((x - 1)^2)/(2 * theta * x)) * (((x - 1)^2) * (2 * x) * 
  (2 * (2 * x * (2 * theta * x)))/((2 * theta * x)^2)^2))/sqrt(2 * 
   pi * theta * (x^3)) - exp(-((x - 1)^2)/(2 * theta * x)) * 
  (((x - 1)^2) * (2 * x)/(2 * theta * x)^2) * (0.5 * (2 * pi * 
  (x^3) * (2 * pi * theta * (x^3))^-0.5))/sqrt(2 * pi * theta * 
   (x^3))^2 - ((exp(-((x - 1)^2)/(2 * theta * x)) * (((x - 1)^2) * 
   (2 * x)/(2 * theta * x)^2) * (0.5 * (2 * pi * (x^3) * (2 * 
  pi * theta * (x^3))^-0.5)) - exp(-((x - 1)^2)/(2 * theta * 
   x)) * (0.5 * (2 * pi * (x^3) * ((2 * pi * theta * (x^3))^-(0.5 + 
  1) * (0.5 * (2 * pi * (x^3)))))))/sqrt(2 * pi * theta * (x^3))^2 - 
   exp(-((x - 1)^2)/(2 * theta * x)) * (0.5 * (2 * pi * (x^3) * 
   (2 * pi * theta * (x^3))^-0.5)) * (2 * (0.5 * (2 * pi * 
  (x^3) * (2 * pi * theta * (x^3))^-0.5) * sqrt(2 * pi * 
  theta * (x^3))))/(sqrt(2 * pi * theta * (x^3))^2)^2)
}

################################################################################
# Log-normal mixture

# 
# A bimodal log-normal mixture, given by the variance and distance between means
# 
dlognormalmix_r <- function(x, theta) {
  V <- theta[1]/2
  M <- theta[2]
  logsd1 <- sqrt( log( (sqrt(4*V + 1) + 1)/2 ) )
  mean1 <- exp(logsd1^2 / 2)
  mean2 <- mean1 + M
  logsd2 <- sqrt( log( V/mean2^2 + 1 ) )
  logmean2 <- log(mean2) - logsd2^2/2
  dlnorm(x, 0, logsd1) + dlnorm(x, logmean2, logsd2)
}

rlognormalmix_r <- function(n, theta) {
  if (theta[1] == 0)
    return(rep(1, n))
  
  V <- theta[1]/2
  M <- theta[2]
  logsd1 <- sqrt( log( (sqrt(4*V + 1) + 1)/2 ) )
  mean1 <- exp(logsd1^2 / 2)
  mean2 <- mean1 + M
  logsd2 <- sqrt( log( V/mean2^2 + 1 ) )
  logmean2 <- log(mean2) - logsd2^2/2
  c(rlnorm(n/2, 0, logsd1) + rlnorm(n/2, logmean2, logsd2))
}

vlognormalmix <- function(theta) {
  theta[1]
}

################################################################################
# Cluster size distributions
################################################################################
# K-truncated Poisson

# Truncated poisson distribution
# Sample from a k-truncated Poisson, where support is strictly greater than k
# modifed from: https://stat.ethz.ch/pipermail/r-help/2005-May/070680.html
rtpois <- function(n, lambda, k=0) {
  qpois(runif(n, sum(dpois(0:k, lambda)), 1), lambda)
}

# Expected value of the k truncated Poisson
etpois <- function(lambda, k=0) {
  one2k <- 1:k
  zero2k <- 0:k
  
  if (k > 0)
    numer_factor2 <- sum(lambda^one2k / factorial(one2k - 1))
  else
    numer_factor2 <- 0
  
  numer <- lambda - exp(-lambda) * numer_factor2
  denom <- 1 - exp(-lambda) * sum(lambda^zero2k / factorial(zero2k))
  numer/denom
}

################################################################################
# Truncated Zeta distribution
# Density, random generator, and expected value of a truncated discrete Pareto 
# (or zeta) distrubition, where support is from xmin to xmax.
# Density
dtzeta <- function(x, alpha, xmin=0, xmax=1e4){
  out <- rep(0, length(x))
  denom <- sum((1:(xmax-xmin))^-alpha / zeta(alpha))
  idx <- (x > xmin)&(x <= xmax)
  out[idx] <- ((x[idx]-xmin)^-alpha / zeta(alpha))/denom
  out
}

# The "easy" way of sampling the truncated Zeta
rtzeta <- function(n, alpha, xmin=0, xmax=1e4){
  sample(x=xmin:xmax, size=n, replace=TRUE,
         prob=dtzeta(x=xmin:xmax, alpha, xmin, xmax))
}

etzeta <- function(alpha, xmin=0, xmax=1e4) {
  k <- 1:(xmax - xmin)
  1/zeta(alpha) * sum(1/(k^(alpha - 1)))
}

################################################################################
# frailty distribution functions can be accessed from these lists

# densities
dfrailty <- list(
  gamma=dgamma_r,
  pvf=dpvf_r,
  posstab=dposstab_r,
  invgauss=dinvgauss_r,
  lognormal=dlognormal_r
)

# Laplace transforms
ltfrailty <- list(
  gamma=lt_dgamma_c,
  pvf=lt_dpvf_c,
  posstab=lt_dposstab_c,
  invgauss=lt_dinvgauss_c,
  lognormal=lt_dlognormal_c
)

# simulation
rfrailty <- list(
  gamma=rgamma_r,
  pvf=rpvf_r,
  posstab=rposstab_r,
  invgauss=rinvgauss_r,
  lognormal=rlognormal_r
)

# variance
vfrailty <- list(
  gamma=vgamma_r,
  pvf=vpvf_r,
  lognormal=vlognormal_r,
  invgauss=vinvgauss_r
)

# initial params for estimation
init.frailty <- list(
  gamma=1,
  pvf=0.5,
  invgauss=1,
  lognormal=1
)

# Soft lower and upper bounds (for estimation)
lb.frailty <- list(
  gamma=1e-2,
  pvf=1e-5,
  posstab=1e-5,
  invgauss=1e-2,
  lognormal=1e-2
)

ub.frailty <- list(
  gamma=Inf,
  pvf=1,
  posstab=1-1e-5,
  invgauss=Inf,
  lognormal=Inf
)

# Hard lower and upper bounds (for generation)
lb.hard.frailty <- list(
  gamma=0,
  pvf=0,
  posstab=0,
  invgauss=0,
  lognormal=0
)

ub.hard.frailty <- list(
  gamma=Inf,
  pvf=1,
  posstab=1,
  invgauss=Inf,
  lognormal=Inf
)