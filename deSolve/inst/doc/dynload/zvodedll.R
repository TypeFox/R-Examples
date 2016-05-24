## =============================================================================
## Implements the test model, as given in the dvode code.
## before trying this code, the FORTRAN program has to be compiled
## this can be done in R:
## system("R CMD SHLIB zvodedll.f")
## do make sure that these files are in the working directory...
## (if not, use setwd() )
## =============================================================================

## the example in "zvode.f",
##
## df/dt = 1i*f    
## dg/dt = -1i*g*g*f 
##
## Initial values are
## g(0) = 1/2.1 and
## z(0) = 1  (same as above)
##
## The analytical solution is
## f(t) = exp(1i*t)   (same as above)
## g(t) = 1/(f(t) + 1.1)

library(deSolve)

## -----------------------------------------------------------------------------
## implementation in R
## -----------------------------------------------------------------------------

ZODE2 <- function(Time, State, Pars) {
  with(as.list(State), {
    df <- 1i * f
    dg <- -1i*g*g*f
    return(list(c(df, dg)))
  })
}

pars  <- NULL
yini  <- c(f = 1+0i, g = 1/2.1+0i)
times <- seq(0, 2*pi, length = 100)

print(system.time(
  out     <- zvode(func = ZODE2, y = yini, parms = NULL, times = times,
               atol = 1e-10, rtol = 1e-10)
))

analytical  <- cbind(f = exp(1i*times), g = 1/(exp(1i*times)+1.1))

#compare numerical solution and the two analytical ones:
tail(cbind(out[, 2], analytical[, 1]))

#----------------------
# the Jacobian:
#----------------------

jac  <- function (t, Y, parameters) {
  PD[2, 2] = -2.0*1i*Y[1]*Y[2]
  PD[2, 1] = -1i*Y[2]*Y[2]
  PD[1, 2] = 0.
  PD[1, 1] = 1i
  return(PD)
} 

print(system.time(
  out2    <- zvode(func = ZODE2, jacfunc = jac, y = yini, parms = NULL,
               times = times, atol = 1e-10, rtol = 1e-10)
))
tail(cbind(out2[, 2], analytical[, 1]))

## -----------------------------------------------------------------------------
## implementation in FORTRAN
## -----------------------------------------------------------------------------

# compiled within R with: system("R CMD SHLIB zvodedll.f")

dyn.load(paste("zvodedll", .Platform$dynlib.ext, sep = ""))
print("FORTRAN DLL passed to zvode")
print(system.time(
  outF <- zvode(func = "fex", jacfunc = "jex", y = yini, parms = NULL,
            times = times, atol = 1e-10, rtol = 1e-10, dllname = "zvodedll",
            initfunc = NULL)
))
tail(cbind(outF[, 2], analytical[, 1]))

