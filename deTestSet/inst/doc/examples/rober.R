## =============================================================================
##
## Rober problem, chemical pyrolysis
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 3
##
## =============================================================================

require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# initial conditions of state variables
yini <- c(1, 0, 0)

# parameters
parms <- c(k1 = 0.04, k2 = 3e7, k3 = 1e4)

# derivative function
Rober <- function(t, y, parms) {
  with (as.list(parms), {
  
      dy1<- -k1*y[1]                + k3*y[2]*y[3]
      dy2<-  k1*y[1] - k2*y[2]*y[2] - k3*y[2]*y[3]
      dy3<-  k2*y[2]*y[2]
      list(c(dy1, dy2, dy3))
  })
}

# -------------------------------------------------------
# run at high resolution 
# -------------------------------------------------------

times <- 10^(seq(from = -5, to = 11, by = 0.1))

# lsoda!
print (system.time(
out <- ode(func = Rober, parms = parms, y = yini,
           times = times, atol = 1e-15, rtol = 1e-15,
           maxsteps = 1e5)
))

print (system.time(
out2 <- ode(func = Rober, parms = parms, y = yini,
           times = times, atol = 1e-15, rtol = 1e-15,
           maxsteps = 1e5, method = mebdfi)
))

print (system.time(
out3 <- ode(func = Rober, parms = parms, y = yini,
           times = times, atol = 1e-15, rtol = 1e-15,
           maxsteps = 1e5, method = gamd)
))


plot(out, out2, out3, lwd = 2, log = "x")
mtext(side = 3, outer = TRUE, line = -1.5, cex = 1.5, "Rober")
