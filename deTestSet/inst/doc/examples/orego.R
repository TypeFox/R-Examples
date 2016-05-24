## =============================================================================
##
## Oregonator problem, chemistry
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
yini <- 1:3

# derivative function
orego <- function(t,y,parms) {
   f1 <- 77.27 * (y[2] + y[1] - y[1]*y[2] - 8.375e-6*y[1]*y[1])
   f2 <- (y[3] - (1. +y[1])*y[2])/77.27
   f3 <- 0.161 *(y[1] - y[3])
   list(c(f1, f2, f3))
}

# -------------------------------------------------------
# run at high resolution 
# -------------------------------------------------------

times <- 0:360
print (system.time(
out <- ode(func = orego, parms = parms, y = yini, times = times,
           atol = 1e-10, rtol = 1e-10, maxsteps = 1e5)
))
print (system.time(
out2 <- ode(func = orego, parms = parms, y = yini, times = times,
           atol = 1e-10, rtol = 1e-10, maxsteps = 1e5, method = mebdfi)
))
print (system.time(
out3 <- ode(func = orego, parms = parms, y = yini, times = times,
           atol = 1e-10, rtol = 1e-10, maxsteps = 1e5, method = gamd)
))
## CRASHES R
print (system.time(
out4 <- ode(func = orego, parms = parms, y = yini, times = times,
            atol = 1e-10, rtol = 1e-10, maxsteps = 1e5, method = cashkarp,
            stiffness = 4)
))

plot(out, out2, out3, lwd = 2, log = "y")
mtext(side = 3, outer = TRUE, line = -1.5, cex = 1.5, "oregonator")

# -------------------------------------------------------
# Compare with exact solution
# -------------------------------------------------------
exact <-c(0.1000814870318523e1, 0.1228178521549917e4, 0.1320554942846706e3)
out[nrow(out),-1]-exact

