#dyn.load("../src/R2Cuba.so");source("../R/divonne.R");source("../R/commoncuba.R")
library("R2Cuba")
# DEMO EXAMPLES
# +++++++++++++++++++++++++++++++
integrand <- function(arg, phase) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
# cat("PHASE", phase)
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # fin integrand

NDIM <- 3
NCOMP <- 1
NMAX <- 4
divonne(NDIM, NCOMP,  integrand, flags=list(verbose=3))
# ------------------------------------------------------------
# Essai avec des arguments ...
integrandp <- function(arg, phase, ...) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
#print("i", ...)
  ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # fin integrandp
divonne(NDIM, NCOMP,  integrandp, quote=FALSE, flags=list(verbose=0))


# ------------------------------------------------------------
# Essai en déplacant les bornes
integrand2 <- function(arg,phase, a) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x-3)*cos(y-2)*exp(z-1);
return(ff)
} # End integrand2
divonne(NDIM, NCOMP,  integrand2, lower=c(3,2,1), upper=c(4,3,2), flags=list(verbose=3))
# ------------------------------------------------------------
# Essai avec mersenne.seed
EPSREL <- 1e-3
EPSABS <- 1e-12
divonne(NDIM, NCOMP, integrand,
  rel.tol=EPSREL,  abs.tol=EPSABS,
      flags=list(verbose=0, pseudo.random=1, mersenne.seed=10))
# ---------------------------------------------------------------
# Essai d'intégrande a pls composantes
integrand3 <- function(arg, phase,a) {
  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
  ff <- sin(x)*cos(y)*exp(z);
  gg <-  1/(3.75 - cos(pi*x) - cos(pi*y) - cos(pi*z));
return(c(ff,gg))
} # End integrand3
divonne(NDIM,2, integrand3,  flags=list(verbose=3))


# ------------------------------------------------------------
# Test de peakfinder: cf peakf.R
# ------------------------------------------------------------
# Essai de key1 différentes
VERBOSE <- 0
KEY1 <- 47
cat(" KEY1=47 pseudo.random=1\n")
divonne(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL, abs.tol=EPSABS,
               flags=list(pseudo.random=1, verbose=VERBOSE) ,
                key1= KEY1)

cat(" KEY1=447 pseudo.random=1\n")
divonne(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL, abs.tol=EPSABS,
               flags=list(pseudo.random=1, verbose=VERBOSE) ,
                key1= 447)


cat(" KEY1=-447 pseudo.random=1\n")
divonne(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL, abs.tol=EPSABS,
               flags=list(pseudo.random=1, verbose=VERBOSE) ,
                key1= -447)


cat(" KEY1=-447 pseudo.random=0\n")
divonne(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
               flags=list(pseudo.random=0, verbose=VERBOSE) ,
                key1= -447)


cat(" KEY1=0 pseudo.random=0\n")
divonne(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
               flags=list(pseudo.random=0, verbose=VERBOSE) ,
                key1= 0)
cat(" KEY1=0,  pseudo.random=1\n")
divonne(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
               flags=list(pseudo.random=1, verbose=VERBOSE) ,
                key1= 0)
