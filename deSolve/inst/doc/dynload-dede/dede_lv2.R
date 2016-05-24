### Lotka-Volterra system with delay

library(deSolve)

derivs <- function(t, y, parms) {
  with(as.list(c(y, parms)), {
    if (t < max(tau1, tau2))
      ytau <- c(1, 1)
    else {
      ytau <- c(
        lagvalue(t - tau1, 1),
        lagvalue(t - tau2, 2)
      )
    }
    dN <- f * N - g * N * P
    dP <- e * g * ytau[1] * ytau[2] - m * P
    list(c(dN, dP), tau=ytau[1], tau=ytau[2])
  })
}

yinit <- c(N=1, P=1)
times <- seq(0, 500)
parms <- c(f=0.1, g=0.2, e=0.1, m=0.1, tau1 = 0.2, tau2 = 50)

## one single run
system.time(
  yout <- dede(y = yinit, times = times, func = derivs, parms = parms)
)
if (!interactive()) pdf(file="dede_lf2.pdf")

plot(yout)

system("R CMD SHLIB dede_lv2.c")
dyn.load(paste("dede_lv2", .Platform$dynlib.ext, sep=""))

## 100 runs
system.time( for (i in 1:100)
  yout2 <- dede(yinit, times = times, func = "derivs", parms = parms,
    dllname = "dede_lv2", initfunc = "initmod", nout = 2)
)

## version "derivs2" (different if tau1 != tau2; respects individual tau
system.time( for (i in 1:100)
  yout3 <- dede(yinit, times = times, func = "derivs2", parms = parms,
    dllname = "dede_lv2", initfunc = "initmod", nout = 2)
)

plot(yout2, yout3) # identical if tau1=tau2


dyn.unload(paste("dede_lv2", .Platform$dynlib.ext, sep=""))

# should be zero
summary(as.vector(yout) - as.vector(yout2))

# can be different from zero
summary(as.vector(yout) - as.vector(yout3))

##
## Fortran Example
##

system("R CMD SHLIB dede_lv2F.f dedeUtils.c")
dyn.load(paste("dede_lv2F", .Platform$dynlib.ext, sep=""))

## 100 runs
system.time( for (i in 1:100)
  yout4 <- dede(yinit, times = times, func = "derivs", parms = parms,
    dllname = "dede_lv2F", initfunc = "initmod", nout = 2)
)

## version "derivs2" (different if tau1 != tau2; respects individual tau
system.time( for (i in 1:100)
  yout5 <- dede(yinit, times = times, func = "derivs2", parms = parms,
    dllname = "dede_lv2F", initfunc = "initmod", nout = 2)
)

plot(yout4, yout5) # identical if tau1=tau2


dyn.unload(paste("dede_lv2F", .Platform$dynlib.ext, sep=""))

# should be zero
summary(as.vector(yout) - as.vector(yout4))

# can be different from zero
summary(as.vector(yout) - as.vector(yout5))

if (!interactive()) dev.off()
