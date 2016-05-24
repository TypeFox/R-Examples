## =============================================================================
## This script demonstrates the performance gain
## from using compiled versions of the ODE models
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package:
##     citation(package="deSolve"       #  doi:10.18637/jss.v033.i09
##     citation(package="FME")          #  doi:10.18637/jss.v033.i03
##     citation(package="growthrates")
## Repository:   https://github.com/tpetzoldt/growthrates
## Mailing list: https://stat.ethz.ch/mailman/listinfo/r-sig-dynamic-models
## =============================================================================

library("growthrates")

## R versions of the models (the C versions are in the package)
grow_twostep.R <- function(time, parms, ...) {
 y0    <- parms[c("yi", "ya")]
 parms <- parms[c("kw", "mumax", "K")]
 out  <-  ode(y0, time, ode_twostep, parms, ...)
}

grow_genlogistic.R <- function(time, parms, ...) {
  y0    <- c(y = unname(parms[c("y0")]))
  parms <- parms[c("mumax", "K", "alpha", "beta", "gamma")]
  out  <-  as.matrix(ode(y0, time, ode_genlogistic, parms, ...))
}


### Two-Step growth model (2 ODE equations)
## R code
system.time(for (i in 1:100)
  o1 <-
    grow_twostep.R(0:100, c(
      yi = 0.01, ya = 0.0, kw = 0.1,	mumax = 0.2, K = 0.1
    )))

## compiled C code
system.time(for (i in 1:100)
  o2 <-
    grow_twostep(0:100, c(
      yi = 0.01, ya = 0.0, kw = 0.1,	mumax = 0.2, K = 0.1
    )))

### extended logistic growth model (2 ODE equations)
## R code
system.time(for (i in 1:100)
  o3 <-
    grow_genlogistic.R(0:100, c(
      y0 = 0.1, mumax = 0.5, K = 10, alpha = 1.2, beta = 1.2, gamma = 1.2
    )))

## compiled C code
system.time(for (i in 1:100)
  o4 <-
    grow_genlogistic(0:100, c(
      y0 = 0.1, mumax = 0.5, K = 10, alpha = 1.2, beta = 1.2, gamma = 1.2
    )))


## check if results are identical
summary(o1 - o2)

summary(o3 - o4)



