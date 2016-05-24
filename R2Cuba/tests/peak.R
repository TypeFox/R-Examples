#dyn.load("../src/cuba.so");source("../R/divonne.R");source("../R/commoncuba.R")
library("R2Cuba")
# Test de peakfinder
# +++++++++++++++++++++++++++++++
NDIM <- 3
NCOMP <- 1
NMAX <- 4
integrand <- function(arg, phase) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
# cat("PHASE", phase)
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # fin integrand

peakf <- function(bounds) {
#  print(bounds)
  x <- matrix(0, ncol=NMAX, nrow=NDIM)
   pas <- 1/(NMAX-1)
   # 1ier point
   x[,1] <- rep(0, NDIM)
   # Les autres points
   for (i in 2:NMAX) {
      x[,i] <- x[,(i-1)] + pas
    }
  return(x)
} #end peakf

print(divonne(NDIM, NCOMP, integrand, flags=list(verbose=0) , peakfinder=peakf, nextra=NMAX))


# -----------------------------------------------------------
# Essai en déplacant les bornes
integrand2 <- function(arg,phase, a) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x-3)*cos(y-2)*exp(z-1);
return(ff)
} # End integrand2

peakf <- function(bounds) {
#  print(bounds)
   x <- matrix(0, ncol=NMAX, nrow=NDIM)
   pas <- 1/(NMAX-1)
   # 1ier point
   x[,1] <- c(3,2,1)
   # Les autres points
   for (i in 2:NMAX) {
      x[,i] <- x[,(i-1)] + pas
    }
   
  return(x)
} #end peakf

print(divonne(NDIM, NCOMP, integrand2, lower=c(3,2,1), upper=c(4,3,2),
               flags=list(verbose=0) ,
                peakfinder=peakf, nextra=NMAX))

  
