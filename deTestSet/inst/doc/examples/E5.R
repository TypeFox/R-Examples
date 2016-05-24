## =============================================================================
##
## E5 problem, chemical pyrolysis      
##
## This code is derived from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##     ODE of dimension 4
##
## =============================================================================

require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# initial conditions of state variables
yini <- c(1.76e-3, rep(1e-20, 3))

# parameters
parms <- c(A = 7.89e-10, B = 1.1e7, C = 1.13e3, M = 1e6)

# derivative function
E3 <- function(t, y, parms) {
  with (as.list(parms), {
  
      dy1<- -A*y[1] - B*y[1]*y[3]
      dy2<-  A*y[1]               - M*C*y[2]*y[3]
      dy3<-  A*y[1] - B*y[1]*y[3] - M*C*y[2]*y[3] + C*y[4]
      dy4<-           B*y[1]*y[3]                 - C*y[4]
      list(c(dy1, dy2, dy3, dy4))
  })
}

# jacobian function
jac <- function(t, y, parms) {   #d( dy/dti)/dyj
  with (as.list(parms), {
     JAC <- matrix (nrow = 4, ncol = 4, 0)
      JAC[1,1]<- -A-B*y[3]
      JAC[1,3]<- -B*y[1]
  
      JAC[2,1] <-  A
      JAC[2,2] <-  -M*C*y[3]
      JAC[2,3] <-  -M*C*y[2]

      JAC[3,1] <-  A - B*y[3]
      JAC[3,2] <-  -M*C*y[3]
      JAC[3,3] <-  -B*y[1] -M*C*y[2]
      JAC[3,4] <-  C

      JAC[4,1] <-  B*y[3]
      JAC[4,3] <-  B*y[1]
      JAC[4,4] <-  - C

      JAC
  })
}

times <- c(0, 10^(seq(-5, 12, by = 0.1)))

print (system.time(
out <- ode(func = E3, jacfunc = jac, parms = parms,
           y = yini, times = times,
           atol = 1e-15, rtol = 1e-15, maxsteps = 1e5, method = "lsode")
))

print (system.time(
out2 <- mebdfi(func = E3, jacfunc = jac, parms = parms,
           y = yini, times = times,
           atol = 1e-16, rtol = 1e-16, maxsteps = 1e6)
))
print (system.time(

out3 <- gamd(func = E3, jacfunc = jac, parms = parms,
           y = yini, times = times,
           atol = 1e-15, rtol = 1e-15, maxsteps = 1e6)
))

plot(out, out2, out3, lwd = 2, log = "xy")
mtext(side = 3, outer = TRUE, line = -1.5, cex = 1.5, "E5 -  pyrolysis")

# -------------------------------------------------------
# Compare with exact solution
# -------------------------------------------------------

exact <-c(y1 = 0.1152903278711829e-290, y2 = 0.8867655517642120e-22,
          y3 = 0.8854814626268838e-22,  y4 = 0.0000000000000000000)
print (out[nrow(out)  , -1] - exact)
print (out2[nrow(out2), -1] - exact)
print (out3[nrow(out3), -1] - exact)
