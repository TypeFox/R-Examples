## =============================================================================
## The following example shows how to create user-defined growth models
## from a system of differential equations in R.
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================

library("growthrates")

## =============================================================================
## create a "growthmodel" with interfaces compatible to package growthrates
## see ?growthmodel for details
## Note:
##   It is essential to use consistent names for parameters and initial values!
## =============================================================================

ode_K_linear <- function (time, init, parms, ...) {
  with(as.list(c(parms, init)), {
    dy <- mumax * y * (1 - y/K)
    dK <- dK
    list(c(dy, dK), log_y=unname(log(y)))
  })
}


grow_K_linear <- function(time, parms, ...) {
  init    <- parms[c("y0", "K")]           # initial values
  names(init) <- c("y", "K")               # the names of state variables
  parms <- parms[c("mumax", "dK")]         # the "real" ODE model parms
  out <- ode(init, time, ode_K_linear, parms)
  out
}

grow_K_linear <- growthmodel(grow_K_linear, pnames=c("y0", "K", "mumax", "deltaK"))

head(grow_K_linear(time=1:10, c(y0=.1, K=1, mumax=0.1, dK = 0.5)))

## =============================================================================
## Fit the model
## =============================================================================

x <- seq(5, 100, 5)

y <- c(0.1, 2.2, 3.1, 1.5, 8.9, 8, 8.4, 9.8, 9.3, 10.6, 12, 13.6,
       13.1, 13.3, 11.6, 14.7, 12.6, 13.9, 16.9, 14.4)



fit <- fit_growthmodel(grow_K_linear,
         p=c(y0=0.1, mumax=0.2, K=10, dK = .1), time=x, y=y)
plot(fit)
summary(fit)
