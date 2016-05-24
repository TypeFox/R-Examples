## =============================================================================
# Example 6 from Shampine et al.
# 
#   This is Example 2 from M. Kubicek et alia, Test examples for comparison 
#   of codes for nonlinear boundary value problems in ordinary differential 
#   equations, B. Childs et al., eds., Codes for Boundary-Value Problems in 
#   Ordinary Differential Equations, Lecture Notes in Computer Science #76, 
#   Springer, New York, 1979.  This example shows how to deal with a singular
#   coefficient arising from reduction of a partial differential equation to
#   an ODE by symmetry.  Also, for the physical parameters considered here,
#   the problem has three solutions.
## =============================================================================

catalyst <- function (x, y, parms) {
 with (as.list(parms),{ 
  dydx <- c(y[2],0)
  tmp  <- f^2 * y[1] * exp(g*b*(1-y[1])/(1+b*(1-y[1])))
  if (x == 0) dydx[2] <- 1/3*tmp
  else        dydx[2] <-  -(2/x)*y[2] + tmp
   return(list(dydx))
  }) 
}   

#  Define the physical parameters for this problem.
parms <- c(
  f = 0.6 ,
  g = 40  ,
  b = 0.2 )

yini <- c(y = NA, dy = 0)
yend <- c(y = 1, dy = NA)
x    <- seq(0, 1, by = 0.05)

## =============================================================================
## three solutions, found with different initial guesses
## =============================================================================
xguess <- c(0,1)
yguess <- matrix(data = 1, nr = 2, nc = 2)

Sol <- bvptwp(func = catalyst, x = x,
              yini = yini, yend = yend, 
              parms = parms, xguess = xguess, yguess = yguess)
plot(Sol, which = "y", type = "l", ylim = c(0, 1))

yguess <- matrix(data = 0.5, nr = 2, nc = 2)
Sol2 <- bvptwp(func = catalyst, x = x,
               yini = yini, yend = yend, 
               parms = parms, xguess = xguess, yguess = yguess)
lines(Sol2[,c("x", "y")])

yguess <- matrix(data = 0.0, nr=2,nc=2)
Sol3 <- bvptwp(func = catalyst, x = x,
               yini = yini, yend = yend, 
               parms = parms, xguess = xguess, yguess = yguess)
lines(Sol3[,c("x", "y")])
