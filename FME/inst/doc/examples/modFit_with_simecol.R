library(simecol)
library(FME)

data(chemostat)

cs1 <- chemostat
solver(cs1) <- "lsoda"

obstime <- seq(0, 20, 2)
yobs <- data.frame(
  X    = c(10, 26, 120, 197, 354, 577, 628, 661, 654, 608, 642),
  S    = c(9.6, 10.2, 9.5, 8.2, 6.4, 4.9, 4.2, 3.8, 2.5, 3.8, 3.9)
)

## Function modCost requires that simulated and observer data both
## contain a "time" column.
## It is already available in ysim (returned by the solver)
## and we add it to yobs too.

yobs <- cbind(time = obstime, yobs)

Cost <- function(p, simObj, obstime, yobs) {
  whichpar <- names(p)
  parms(simObj)[whichpar] <- p
  times(simObj) <- obstime
  ysim <- out(sim(simObj))
  modCost(ysim, yobs, weight="std")
}

## only methods that converged without constraints in reasonable time ...

system.time({
Fit <- modFit(p = c(vm = 10, km = 10), f = Cost, simObj = cs1,
  obstime = obstime, yobs = yobs, method = "Nelder")
  print(summary(Fit), cov = FALSE)
  deviance(Fit)
})


system.time({
  Fit <- modFit(p = c(vm = 10, km = 10), f = Cost, simObj = cs1,
    obstime = obstime, yobs = yobs, method = "BFGS")
  print(summary(Fit))
  deviance(Fit)
})

system.time({
  Fit <- modFit(p = c(vm = 10, km = 10), f = Cost, simObj = cs1,
    obstime = obstime, yobs = yobs, method = "CG")
  print(summary(Fit))
  deviance(Fit)
})

system.time({
  Fit <- modFit(p = c(vm = 10, km = 10), f = Cost, simObj = cs1,
    obstime = obstime, yobs = yobs, method = "Newton")
  print(summary(Fit))
  deviance(Fit)
})

## the parameter fitting function of simecol for comparison
## it is less powerful and only a little bit faster
system.time({
  whichpar <- c("vm", "km")
  parms(cs1)[whichpar] <- c(vm=5, km=10)
  parms(cs1)
  res <- fitOdeModel(cs1, whichpar = whichpar, obstime = obstime,
                     yobs = yobs[-1], method = "Nelder")
  res$par
})
