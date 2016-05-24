## =============================================================================
## The following example shows how to use compiled growth models
## from inline code, by using the 'cOde' package of Daniel Kaschek
## Note: This example needs the R development tools.
##  - suitable compilers on Linux and Mac
##  - Rtools on Windows from https://cran.r-project.org/bin/windows/Rtools/
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================


library("growthrates")
library("cOde")

## =============================================================================
## define a system of ODEs and compile it
## =============================================================================
ode_K_linear <- funC(c(
  y = "mumax * y * (1-y/K)",
  K = "dK"
))

yini <- c(y = 1, K = 10)
parms = c(mumax = 0.1, dK = 0.05)

## run the model
out1 <- odeC(yini, times=0:100, ode_K_linear, parms = parms)

## generate artificial test data with normal distributed noise
x <- seq(5, 100, 5)
y <- odeC(yini, x, ode_K_linear, parms)[, "y"] + rnorm(x)

## =============================================================================
## create a "growthmodel" with interfaces compatible to package growthrates
## see ?growthmodel for details
## Note:
##   It is essential to use consistent names for parameters and initial values!
## =============================================================================

grow_K_linear <- function(time, parms, ...) {
  init    <- parms[c("y0", "K")]  # initial values
  names(init) <- c("y", "K")      # force names
  out <- odeC(init, time, ode_K_linear, parms)
  cbind(out, log_y = log(out[,"y"]))
}

## convert this to an object, (maybe needed by future extensions)
# grow_K_linear <- growthmodel(grow_K_linear, pnames=c("y0", "mumax", "K", "dK"))

## Test the growthmodel.
## Columns with names 'time', 'y' and 'log_y' are mandatory.
head(grow_K_linear(time=x, c(y0=1, mumax=0.1, K=10, dK = 0.1)))

## =============================================================================
## Fit the model
## =============================================================================
fit <- fit_growthmodel(
  grow_K_linear, p=c(y0=1, mumax=0.1, K=10, dK = 0.1), time=x, y=y)

plot(fit)
summary(fit)

## unload DLL and cleanup, should ideally work in a temp dir
dll <- paste(ode_K_linear, .Platform$dynlib.ext, sep="")
dyn.unload(dll)
unlink(dll)
unlink(paste(ode_K_linear, ".c", sep=""))
unlink(paste(ode_K_linear, ".o", sep=""))






