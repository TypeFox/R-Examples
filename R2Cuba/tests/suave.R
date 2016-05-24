#dyn.load("../src/R2Cuba.so");source("../R/suave.R"); source("../R/commoncuba.R")
library("R2Cuba")
integrand <- function(arg, weight) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
 #  print( weight)
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # end  integrand

NDIM <- 3
NCOMP <- 1
EPSREL <- 1e-3
EPSABS <- 1e-12
VERBOSE <- 2
LAST <- 1

suave(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
      flags=list(verbose=0))
# ----------------------------------------
# Essai avec ...
integrandp <- function(arg, weight, ...) {

  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
 # print( weight, ...)
    ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # fin integrandp
suave(NDIM, NCOMP, integrandp, digits=1,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=0))
# ----------------------------------------

# Essai en déplacant les bornes
integrand2 <- function(arg, weight) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
   ff <- sin(x-3)*cos(y-2)*exp(z-1);
return(ff)
} # End integrand2
suave(NDIM, NCOMP,  integrand2, lower=c(3,2,1), upper=c(4,3,2), flags=list(verbose=3))
