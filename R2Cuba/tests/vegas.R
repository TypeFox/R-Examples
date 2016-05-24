#dyn.load("../src/R2Cuba.so");source("../R/vegas.R"); source("../R/commoncuba.R")
library("R2Cuba")

# DEMO EXAMPLES
# +++++++++++++++++++++++++++++++
integrand <- function(arg, weight) {

  x <- arg[1]
  y <- arg[2]
  z <- arg[3]
#  print( weight)
    ff <- sin(x)*cos(y)*exp(z);
return(ff)
} # fin integrand
# ----------------------------------------

NDIM <- 3
NCOMP <- 1
EPSREL <- 1e-3
EPSABS <- 1e-12
VERBOSE <- 0

vegas(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=VERBOSE))
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
vegas(NDIM, NCOMP, integrandp, digits=1,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=0))
# ----------------------------------------
# Essai avec state
vegas(NDIM, NCOMP, integrand,
                  rel.tol=1e-50,  abs.tol=1e-50,
             flags=list(verbose=VERBOSE), state.file="toto")

vegas(NDIM, NCOMP, integrand,
                   rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=VERBOSE), state.file="toto")

# ----------------------------------------
# Essai avec mersenne
vegas(NDIM, NCOMP, integrand,
  rel.tol=EPSREL,  abs.tol=EPSABS,
      flags=list(verbose=0, pseudo.random=1, mersenne.seed=10))

# ----------------------------------------
# Essai avec nbatch
vegas(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=0, vegas.nbatch=500))
# ----------------------------------------
# Essai avec gridno
#dyn.load("../src/cuba.so");source("../R/vegas.R"); source("../R/commoncuba.R")
vegas(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=3, pseudo.random=1, mersenne.seed=10,
               vegas.gridno=4))
vegas(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=0, pseudo.random=1, mersenne.seed=10,
               vegas.gridno=4))
vegas(NDIM, NCOMP, integrand,
                  rel.tol=EPSREL,  abs.tol=EPSABS,
             flags=list(verbose=0, pseudo.random=1, mersenne.seed=10,
               vegas.gridno=-4))
