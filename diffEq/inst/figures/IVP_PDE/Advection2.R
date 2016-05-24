## =============================================================================
## Flux limiters
## Figure 9.4 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================

adv.func <- function(t, y, p, adv.method)
  list(advection.1D(C = y, C.up = y[N], C.down = y[1], 
                   v = 0.1, adv.method = adv.method, 
                   dx = xgrid)$dC)

v  <- 0.1
xgrid <- setup.grid.1D(0.3, 1.3, N = 50)
x     <- xgrid$x.mid
N     <- length(x)

times <- seq(0, 20, 0.01)
ii    <- 401
T     <- times[ii]
shift <- T*v
yana  <- sin(pi*(x-T*v))^50

yini  <- sin(pi * x)^50

out1 <- ode.1D(y = yini, func = adv.func, times = times, 
              parms = NULL, method = "euler", dimens = N, 
              adv.method = "muscl")
out2 <- ode.1D(y = yini, func = adv.func, times = times, 
              parms = NULL, method = "euler", dimens = N, 
              adv.method = "super")

windows(height = 6, width = 8)

par(mfrow = c(1, 2))
plot(x, out1[1, -1], type = "l", xlab = "x", ylab = "y", 
     main = "muscl", ylim  = c(0, 1))
lines(x, out1[ii, -1], lwd = 2)
lines(x, yana, col = "darkgrey", lty = 1, lwd = 2)
#0.25*v/dx
writelabel("A")

plot(x, out2[1, -1], type = "l", xlab = "x", ylab = "y", 
     main = "superbee", ylim  = c(0, 1))
lines(x, out2[ii, -1], lwd = 2)
lines(x, yana, col = "darkgrey", lty = 1, lwd = 2)
#0.25*v/dx
writelabel("B")
