lookbackTime <-
function(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71){
  ## Function to calculate lookback time from redshift z
  ## Refs:
    # "Principles of Physical Cosmology," P.J. Peebles, Princeton c. 1993, 5.63
    # "First-year WMAP observations...", Spergel et al., ApJS 148:175 (2003)
  ## A. Harris 8/2/2008
  ## updated to handle vector of input values z 5/15/2011
  ## Takes:
  ##   Density parameters: omega.m, omega.lambda
  ##   Hubble constant: H.0 in km/s/Mpc
  ## Returns:
  ##   Lookback time in Gyr
  ## Notes:
  ##   Default parameters from WMAP concordance cosmology
  ##   Inverse problem, age of Earth example:
  ##     uniroot(function(x) lookbackTime(x) - 4.6, c(0,2))$root


  ## Calculate curvature parameter, scaling
  omega.k <- 1 - omega.m - omega.lambda
  scale <- 977.814/H.0  # includes unit conversions from km/s/Mpc to Gyr
  ## Function to integrate
  f <- function(x){
         1./((1+x) * sqrt(omega.m*(1+x)^3 + omega.k*(1+x)^2 + omega.lambda))
         }
  ## Lookback time
  len <- length(z)
  t <- numeric(len)
  for (i in 1:len) t[i] <- scale * integrate(f, 0, z[i])$value
  return(t)
}

