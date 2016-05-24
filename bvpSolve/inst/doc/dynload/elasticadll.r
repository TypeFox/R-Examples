## #############################################################################
## Implements the elastica model as a dll
## before trying this code, the fortran code has to be compiled
## this can be done in R:
## system("R CMD SHLIB elastica.f")
## or system("R CMD SHLIB elasticaC.c") or
## system("gfortran -shared -o elastica.dll elastica.f")
## do make sure that the file is in the working directory...
## (if not, use setwd() )
## #############################################################################

## =============================================================================
## R-functions to solve the Elastica equation
## =============================================================================

#  The F-function, the differential system :

Elastica <- function (x, y, pars) {

  list( c(cos(y[3]),
          sin(y[3]),
          y[4],
          y[5]*cos(y[3]),
          0))
}

# The analytic Jacobian for the F-function:

Jac <- matrix(nr=5,nc=5,0)
Jac[3,4]=1.0
Jac[4,4]=1.0

dfsub <- function (x, y, pars) {
      Jac[1,3]=-sin(y[3])
      Jac[2,3]=cos(y[3])
      Jac[4,3]=-y[5]*sin(y[3])
      Jac[4,5]=Jac[2,3]
      Jac
}

require(bvpSolve)

# The G-function

bound <- function (i, y, pars)  {
    if (i <=2) return(y[i])
    else if (i == 3) return(y[4])
    else if (i == 4) return(y[2])
    else if (i == 5) return(y[3]+pi/2)
}

# The analytic Jacobian for the G-function:
jacbound <- function(i, y, pars)  {
    JJ <- rep(0,5)
         if (i <= 2) JJ[i] =1.0
    else if (i == 3) JJ[4] =1.0
    else if (i == 4) JJ[2] =1.0
    else if (i == 5) JJ[3] =1.0
    JJ
}

## =============================================================================
## Different ways to solve the same model
## =============================================================================

niter <- 100
# Only the F-function is supplied, jacobians estimated by perturbation
print(system.time(
for ( i in 1:niter)
Sol <- bvptwp(func = Elastica,
              yini = c(x = 0,  y = 0, p = NA,   k = 0,  F = NA),
              yend = c(x = NA, y = 0, p = -pi/2,k = NA, F = NA),
              x = seq(0, 0.5, len = 16) )
)/niter)

# the model with the F-function and analytical F-jacobian

print(system.time(
for ( i in 1:niter)
Sol2 <- bvptwp(func=Elastica, jacfunc = dfsub,
              yini = c(x = 0, y = 0, p = NA,    k = 0,  F = NA),
              yend = c(x = NA,y = 0, p = -pi/2, k = NA, F = NA),
              x = seq(0, 0.5, len = 16) )
)/niter)

# the model with the F-function, analytical F- and G-jacobian
print(system.time(
for ( i in 1:niter)
Sol3 <- bvptwp(func = Elastica, jacfunc = dfsub, jacbound = jacbound,
              yini = c(x = 0,  y = 0, p = NA,    k = 0,  F = NA),
              yend = c(x = NA, y = 0, p = -pi/2, k = NA, F = NA),
              x = seq(0, 0.5, len = 16) )
)/niter)

# the model with the F- and G-function, analytical F- and G-jacobian
print(system.time(
for ( i in 1:niter)
Sol4 <- bvptwp(leftbc = 3, ncomp = 5,
              func = Elastica, jacfunc = dfsub, bound = bound, jacbound = jacbound,
              x = seq(0, 0.5, len = 16) )
)/niter)

# This model cannot be solved with the shooting method...

# the model implemented in FORTRAN
dyn.load("elastica.dll")

print(system.time(
for (i in 1:niter)
outF <- bvpcol(ncomp = 5,
               x = seq(0, 0.5, len = 16), leftbc = 3,  
               func = "fsub", jacfunc = "dfsub", bound = "gsub",
               jacbound = "dgsub", xguess = Sol4[,1], yguess = t(Sol4[,-1]),
               dllname = "elastica")
)/niter)

dyn.unload("elastica.dll")

# the model implemented in C
dyn.load("elasticaC.dll")

print(system.time(
for (i in 1:niter)
outC <- bvptwp(ncomp = 5,
               x = seq(0, 0.5, len = 16), leftbc = 3,
               func = "fsub", jacfunc = "dfsub",
               bound = "gsub", jacbound = "dgsub",
               dllname = "elasticaC")
)/niter)

dyn.unload("elasticaC.dll")
# differences with respect to the Fortran solution:

max(abs((outC-outF)))
max(abs((Sol-outF)))
max(abs((Sol2-outF)))
max(abs((Sol3-outF)))
max(abs((Sol4-outF)))


