D.M <-
function(z, omega.m = 0.27, omega.lambda = 0.73, H.0 = 71){
  ## Function to calculate transverse (and line-of-sight) comoving distance
  ## from Hogg (16), also (15) astro-ph/9905116
  ## A. Harris 4/17/2011, 10/3/2011
  ## Takes:
  ##   Density parameters: omega.m, omega.lambda
  ##   Hubble constant: H.0 in km/s/Mpc
  ## Returns:
  ##   Comoving transverse distance (Hogg (16))

  ## Calculate curvature parameter, scaling
  omega.k <- 1 - omega.m - omega.lambda
  #D.H <- 1               # no scaling
  #D.H <- 9.25e25/H.0     # meters
  D.H <- 3.e5/H.0         # Mpc

  ## Function to integrate
  f <- function(x){
         1./(sqrt(omega.m*(1+x)^3 + omega.k*(1+x)^2 + omega.lambda))
         }
  ## D.C
  D.C <- numeric(length(z))
  for (i in 1:length(z)) {
      D.C[i] <- ifelse(is.na(z[i]), NA, D.H * integrate(f, 0, z[i])$value)
  }

  ## D.M
  D.M <- switch(sign(round(omega.k, 4))+2,
            {sqk <- sqrt(omega.k)
             D.H/sqk*sin(sqk*D.C/D.H)},
            {D.C},
            {sinh <- function(x) 0.5*(exp(x)-exp(-x))
             sqk <- sqrt(omega.k)
             D.H/sqk*sinh(sqk*D.C/D.H)})
  return(D.M)
}

