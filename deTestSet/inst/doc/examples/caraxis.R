## =============================================================================
##
##     This file is adapted from the Test Set for IVP solvers
##     http://www.dm.uniba.it/~testset/
##
##        Car Axis problem (in index 3 formulation)
##        index 3 DAE of dimension 10
##
##     This is revision
##     $Id: caraxis.F,v 1.2 2006/10/02 10:29:13 testset Exp $
##
## =============================================================================
require(deTestSet)

# -------------------------------------------------------
# problem formulation
# -------------------------------------------------------

# residual function
caraxis <- function(t, y, pars){
  with(as.list(c(pars,y)), {
      f <- rep(0, 10)

      yb  <- r*sin(w*t)
      xb  <- sqrt(L*L - yb*yb)
      Ll  <- sqrt(xl^2 + yl^2)
      Lr  <- sqrt((xr-xb)^2 + (yr-yb)^2)

      f[1:4] <- y[5:8]
      k <- M*eps*eps/2

      f[5]  <- (L0-Ll)*xl/Ll  + lam1*xb + 2*lam2*(xl-xr)
      f[6]  <- (L0-Ll)*yl/Ll  + lam1*yb + 2*lam2*(yl-yr) - k*g
      f[7]  <- (L0-Lr)*(xr-xb)/Lr - 2*lam2*(xl-xr)
      f[8]  <- (L0-Lr)*(yr-yb)/Lr - 2*lam2*(yl-yr)-k*g

      f[9]  <- xb*xl+yb*yl
      f[10] <- (xl-xr)^2+(yl-yr)^2-L*L
      
      list( f )
  })
}


# parameters
pars <- c(eps = 1e-2, M = 10, L = 1, L0 = 0.5,
          r   = 0.1,  w = 10, g = 1)

k    <- with (as.list(pars), M*eps*eps/2)
mass <- diag (nrow = 10, c(rep(1, 4), rep (k, 4), 0, 0))

# initial conditions: state variables
yini <-  with (as.list(pars),
   c(xl  = 0, yl  = L0,    xr  = L, yr = L0, xla = -L0/L,
     yla = 0, xra = -L0/L, yra = 0, lam1 = 0, lam2 = 0)
              )

# initial conditions: derivates
dyini <- rep(0, times = 10)
FF    <- caraxis(0, yini, pars)[[1]]
dyini[1:4] <- yini[5:8]
dyini[5:8] <- 2/pars["M"]/(pars["eps"])^2 * FF[5:8]

# check consistency of initial condition: delt should be = 0.
mass %*% dyini - caraxis(0, yini, pars)[[1]]

# -------------------------------------------------------
# run at high resolution 
# -------------------------------------------------------
times <- seq(0, 3, by = 0.01)
nind  <- c(4, 4, 2)   # index 1, 2 and 3 variables
print (system.time(
out   <- mebdfi(y = yini, dy = dyini, times, func = caraxis,
             mass = mass, parms = pars, nind = nind,
             rtol = 1e-10, atol = 1e-10)
))
print (system.time(
out2   <- gamd(y = yini, times, func = caraxis,
             mass = mass, parms = pars, nind = nind,
             rtol = 1e-10, atol = 1e-10)
))
print (system.time(
out3   <- radau(y = yini, times, func = caraxis,
             mass = mass, parms = pars, nind = nind,
             rtol = 1e-10, atol = 1e-10)
))

plot(out, out2, out3, which = 1:4, lwd = 2, ask = FALSE)

diagnostics(out)
