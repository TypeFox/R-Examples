## =============================================================================
## The advection equation
## Figure 9.3 from Soetaert, Cash and Mazzia, 
## Solving differential equations in R
## =============================================================================
require(ReacTran)

upwind <- function(t, y, p)
  list(tran.1D(C = y, C.up = y[N],  C.down = y[1], 
               v = v, AFDW = 1, dx=dx)$dC)

forward <- function(t, y, p)
  list(tran.1D(C = y, C.up = y[N], C.down = y[1], 
               v = v, AFDW = 0, dx = dx)$dC)

central <-function(t, y, p)
  list(tran.1D(C = y, C.up = y[N], C.down = y[1], 
               v = v, AFDW = 0.5, dx = dx)$dC)

limit <-function(t, y, p, met = "p3")
  list(advection.1D(C = y, C.up = y[N], C.down = y[1], 
                    v = v, adv.method = met, dx = dx)$dC)

v  <- 0.1
dx <- 1/50
x  <- seq(0.3, 1.3, by = dx)
N  <- length(x)

y0   <- 0
yini <-  sin(pi*x)^50
times <- seq(0, 20, 0.01)
ii    <- 401
T     <- times[ii]
shift <- T*v
yana  <- sin(pi*(x-T*v))^50

# 0.01 * v/dx
par(mfrow = c(3, 2))
out1 <- ode.1D(y = yini, func = upwind, times = times, 
               parms = NULL, method = "euler", dimens = N)
plot(x, out1[1,-1], type = "l", xlab = "x", ylab = "y", 
     main = "1st order upwind, c = 0.05")
lines(x, out1[ii,-1], lwd = 2)
lines(x, yana, col = "darkgrey", lty = 2, lwd = 2)
writelabel("A")

out3 <- ode.1D(y = yini, func = limit, times = times, 
               parms = NULL, method = "euler", dimens = N)
plot(x, out3[1,-1], type = "l", xlab="x", ylab = "y", 
     main = "3-rd order upwind")
lines(x, out3[ii, -1], lwd = 2)
lines(x, yana, col = "darkgrey", lty = 2, lwd = 2)
writelabel("B")

out2 <- ode.1D(y = yini, func = central, times = times, 
               parms = NULL, method = "euler", dimens = N)
plot(x, out2[1, -1], type = "l", ylim = c(-0.4,1), xlab = "x", 
     ylab = "y", main = "centered")
lines(x, out2[ii, -1], lwd = 2)
lines(x, yana, col = "darkgrey", lty = 2, lwd = 2)
writelabel("C")

plot(1, type = "n", xlab="", ylab = "", axes = FALSE)
legend("center", lty = c(1,1,2), col = c("black", "black", "grey"), 
       lwd = c(1, 2, 2),
  legend=c("initial", "numerical at t=4", "analytical at t=4"))


dx <- 1/100
x <- seq(0.3, 1.3, by = dx)
N <- length(x)

y0 <- 0
yini <-  sin(pi*x)^50
Times <- seq(0, 20, 0.1)
ii <- which(Times == T)

shift <- T*v
yana <- sin(pi*(x-T*v))^50
out1 <- ode.1D(y = yini, func = upwind, times = Times, 
               parms = NULL, method = "euler", dimens = N)
plot(x, out1[1, -1], type = "l", xlab = "x", ylab = "y", 
     main = "1st order upwind, c = 1")
lines(x, out1[ii, -1], lwd = 2)
lines(x, yana , col = "darkgrey", lty = 2, lwd = 2)
writelabel("D")


dx <- 1/50
x <- seq(0.3, 1.3, by = dx)
N <- length(x)

y0 <- 0
yini <-  sin(pi*x)^50
Times <- seq(0, 20, 0.25)
ii <- which(Times == T)

shift <- T*v
yana <- sin(pi*(x-T*v))^50
out1 <- ode.1D(y = yini, func = upwind, times = Times, 
               parms = NULL, method = "euler", dimens = N)
plot(x, out1[1, -1], type = "l", xlab = "x", ylab = "y", 
     main = "1st order upwind, c = 1.25", ylim  = c(-1, 2))
lines(x, out1[ii, -1], lwd = 2)
lines(x, yana, col = "darkgrey", lty = 2, lwd = 2)
#0.25*v/dx
writelabel("E")

