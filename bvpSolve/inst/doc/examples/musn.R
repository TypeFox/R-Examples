## =============================================================================
##   This is the example for MUSN in U. Ascher, R. Mattheij, and R. Russell,
##   Numerical Solution of Boundary Value Problems for Ordinary Differential
##   Equations, SIAM, Philadelphia, PA, 1995.  MUSN is a multiple shooting
##   code for nonlinear BVPs.  The problem is
##
##      u' =  0.5*u*(w - u)/v
##      v' = -0.5*(w - u)
##      w' = (0.9 - 1000*(w - y) - 0.5*w*(w - u))/z
##      z' =  0.5*(w - u)
##      y' = -100*(y - w)
##
##   The interval is [0 1] and the boundary conditions are
##
##      u(0) = v(0) = w(0) = 1,  z(0) = -10,  w(1) = y(1)
##
## note: there are two solutions...
## =============================================================================

require(bvpSolve)

## =============================================================================
## First method: shooting
## =============================================================================

# Derivatives

musn <- function(t,Y,pars)  {
  with (as.list(Y),
  {
   du <- 0.5*u*(w-u)/v
   dv <- -0.5*(w-u)
   dw <- (0.9-1000*(w-y)-0.5*w*(w-u))/z
   dz <- 0.5*(w-u)
   dy <- -100*(y-w)
   return(list(c(du, dv, dw, dz, dy)))
  })
}

x<- seq(0,1,by=0.05)

# Residual function for yend...
res  <- function (Y,yini,pars)  with (as.list(Y), w-y)

# Initial values; NA= not available
init <- c(u=1, v=1, w=1, z=-10, y=NA)

print(system.time(
sol   <-bvpshoot(func = musn, yini= init, yend = res, x = x, 
           guess = 1, atol = 1e-10, rtol = 0)
))

# second solution...
sol2  <-bvpshoot(func = musn, yini = init, yend = res, x = x,
           guess = 0.9, atol = 1e-10, rtol = 0)

pairs(sol)

# check the solution by simple integration...
yini <- sol[1,-1]
out  <- as.data.frame(ode(fun = musn, y = yini, parms = 0,
     times = x, atol = 1e-10, rtol = 0))
out$w[nrow(out)]-out$y[nrow(out)]

## =============================================================================
## Solution method 2 : bvptwp 
## =============================================================================
# Does not work unless good initial conditions are used

bound <- function(i, y, pars) {
  with (as.list(y), {
    if (i == 1) return (u-1)
    if (i == 2) return (v-1)
    if (i == 3) return (w-1)
    if (i == 4) return (z+10)
    if (i == 5) return (w-y)
 })
}

xguess <- seq(0, 1, len = 5)
yguess <- matrix(nc = 5,data = rep(c(1,1,1,-10,0.91),5))
rownames(yguess) <- c("u", "v", "w", "z", "y")

print(system.time(
Sol <- bvptwp(yini = NULL, x = x, func = musn, bound = bound,
              xguess = xguess, yguess = yguess, leftbc = 4,
              atol = 1e-10)
))
plot(Sol)


# same using bvpshoot - not so quick
print(system.time(
Sol2 <- bvpshoot(yini = NULL, x = x, func = musn, bound = bound,
       leftbc = 4,  guess=c(u=1, v=1, w=1, z=-10, y=0), atol = 1e-10)
)) 
