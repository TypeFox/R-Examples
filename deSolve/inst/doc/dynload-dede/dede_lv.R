### Simple DDE, adapted version of ?dede example from package deSolve

library(deSolve)

derivs <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    if (t < tau)
      ytau <- c(1, 1)
    else
      ytau <- lagvalue(t - tau, c(1, 2))

    dN <- f * N - g * N * P
    dP <- e * g * ytau[1] * ytau[2] - m * P
    list(c(dN, dP), tau=ytau[1], tau=ytau[2])
  })
}

yinit <- c(N=1, P=1)
times <- seq(0, 500)
parms <- c(f=0.1, g=0.2, e=0.1, m=0.1, tau = .2)

## one single run
system.time(
  yout <- dede(y = yinit, times = times, func = derivs, parms = parms)
)

if(!interactive()) pdf(file="dede_lv.pdf")

plot(yout)


system("R CMD SHLIB dede_lv.c")
dyn.load(paste("dede_lv", .Platform$dynlib.ext, sep=""))

## 100 runs
system.time( for (i in 1:100)
  yout2 <- dede(yinit, times = times, func = "derivs", parms = parms,
    dllname = "dede_lv", initfunc = "initmod", nout = 2)
)

dyn.unload(paste("dede_lv", .Platform$dynlib.ext, sep=""))

plot(yout2, main=c("y", "ytau"))

## Fortran models still need the c code in dedeUtils.c.
## However, as long as you just use the lagvalue() and lagderiv()
## supplied with deSolve, dedeUtils.c works as is.

system("R CMD SHLIB dede_lvF.f dedeUtils.c")
dyn.load(paste("dede_lvF", .Platform$dynlib.ext, sep=""))

## 100 runs
system.time( for (i in 1:100)
  yout3 <- dede(yinit, times = times, func = "derivs", parms = parms,
    dllname = "dede_lvF", initfunc = "initmod", nout = 2)
)

dyn.unload(paste("dede_lvF", .Platform$dynlib.ext, sep=""))

plot(yout3, main=c("y", "ytau"))

if(!interactive()) dev.off()
