### Simple DDE, adapted version of ?dede example from package deSolve

library(deSolve)

derivs <- function(t, y, parms) {
  with(as.list(parms), {
    if (t < tau)
      ytau <- 1
    else
      ytau <- lagvalue(t - tau)

    dy <- k * ytau
    list(c(dy), ytau=ytau)
  })
}

yinit <- c(y=1)
times <- seq(0, 30, 0.1)
parms <- c(tau = 1, k = -1)

## one single run
system.time(
  yout <- dede(y = yinit, times = times, func = derivs, parms = parms)
)

if (!interactive()) pdf(file="dedesimple.pdf")
plot(yout, main = c("dy/dt = -y(t-1)", "ytau"))


system("R CMD SHLIB dedesimple.c")
#dyn.load("dedesimple.dll")
dyn.load(paste("dedesimple", .Platform$dynlib.ext, sep=""))

## 100 runs
system.time( for (i in 1:100)
  yout2 <- dede(yinit, times = times, func = "derivs", parms = parms,
    dllname = "dedesimple", initfunc = "initmod", nout = 1)
)

#dyn.unload("dedesimple.dll")
dyn.unload(paste("dedesimple", .Platform$dynlib.ext, sep=""))

plot(yout2, main=c("y", "ytau"))

## Fortran example

system("R CMD SHLIB dedesimpleF.f dedeUtils.c")
dyn.load(paste("dedesimpleF", .Platform$dynlib.ext, sep=""))

## 100 runs
system.time( for (i in 1:100)
  yout3 <- dede(yinit, times = times, func = "derivs", parms = parms,
    dllname = "dedesimpleF", initfunc = "initmod", nout = 1)
)

#dyn.unload("dedesimple.dll")
dyn.unload(paste("dedesimpleF", .Platform$dynlib.ext, sep=""))

plot(yout3, main=c("y", "ytau"))

if (!interactive()) dev.off()
