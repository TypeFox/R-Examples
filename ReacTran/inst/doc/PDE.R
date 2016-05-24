### R code from vignette source 'PDE.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library(ReacTran)
options(prompt = " ")


###################################################
### code chunk number 2: PDE.Rnw:175-176
###################################################
Grid <- setup.grid.1D(N = 1000, L = 10)


###################################################
### code chunk number 3: PDE.Rnw:182-184
###################################################
r  <- setup.prop.1D(grid = Grid, func = function(r) r)
r2 <- setup.prop.1D(grid = Grid, func = function(r) r^2)


###################################################
### code chunk number 4: PDE.Rnw:194-199
###################################################
library(ReacTran)
pde1D <-function(t, C, parms, A = 1)  {
  tran   <- tran.1D (C = C, A = A, D = D, C.down = Cext, dx = Grid)$dC
  list(tran - Q)  # the return value: rate of change
}


###################################################
### code chunk number 5: PDE.Rnw:204-207
###################################################
D    <- 1    # diffusion constant
Q    <- 1    # uptake rate
Cext <- 20


###################################################
### code chunk number 6: PDE.Rnw:232-237
###################################################
library(rootSolve)
Cartesian   <- steady.1D(y = runif(Grid$N),
  func = pde1D, parms = NULL, nspec = 1, A = 1)
Cylindrical <- steady.1D(y = runif(Grid$N),
  func = pde1D, parms = NULL, nspec = 1, A = r)


###################################################
### code chunk number 7: PDE.Rnw:239-243
###################################################
print(system.time(
  Spherical   <- steady.1D(y = runif(Grid$N),
    func = pde1D, parms = NULL, nspec = 1, A = r2)
))


###################################################
### code chunk number 8: pde
###################################################
plot(Cartesian, Cylindrical, Spherical, grid = Grid$x.mid, 
   main = "steady-state PDE", xlab = "x", ylab = "C",
   col = c("darkgreen", "blue", "red"), lwd = 3, lty = 1:3)

legend("bottomright", c("cartesian", "cylindrical", "spherical"),
  col = c("darkgreen", "blue", "red"), lwd = 3, lty = 1:3)


###################################################
### code chunk number 9: figpde
###################################################
plot(Cartesian, Cylindrical, Spherical, grid = Grid$x.mid, 
   main = "steady-state PDE", xlab = "x", ylab = "C",
   col = c("darkgreen", "blue", "red"), lwd = 3, lty = 1:3)

legend("bottomright", c("cartesian", "cylindrical", "spherical"),
  col = c("darkgreen", "blue", "red"), lwd = 3, lty = 1:3)


###################################################
### code chunk number 10: PDE.Rnw:275-278
###################################################
max(abs(Q/6/D*(r2$mid - 10^2) + Cext - Spherical$y))
max(abs(Q/4/D*(r2$mid - 10^2) + Cext - Cylindrical$y))
max(abs(Q/2/D*(r2$mid - 10^2) + Cext - Cartesian$y))


###################################################
### code chunk number 11: PDE.Rnw:285-291
###################################################
require(deSolve)
times <- seq(0, 100, by = 1)
system.time(
  out <- ode.1D(y = rep(1, Grid$N), times = times, func = pde1D,
    parms = NULL, nspec = 1, A = r2)
)


###################################################
### code chunk number 12: PDE.Rnw:298-299
###################################################
tail(out[, 1:4], n = 3)


###################################################
### code chunk number 13: yy
###################################################
image(out, grid = Grid$x.mid, xlab = "time, days", 
     ylab = "Distance, cm", main = "PDE", add.contour = TRUE)


###################################################
### code chunk number 14: PDE.Rnw:356-360
###################################################
dx    <- 0.2
xgrid <- setup.grid.1D(-100, 100, dx.1 = dx)
x     <- xgrid$x.mid
N     <- xgrid$N


###################################################
### code chunk number 15: PDE.Rnw:363-367
###################################################
uini  <- exp(-0.05 * x^2)
vini  <- rep(0, N)
yini  <- c(uini, vini)
times <- seq (from = 0, to = 50, by = 1)


###################################################
### code chunk number 16: PDE.Rnw:374-382
###################################################
wave <- function (t, y, parms) {
  u1 <- y[1:N]
  u2 <- y[-(1:N)]

  du1 <- u2
  du2 <- tran.1D(C = u1, C.up = 0, C.down = 0, D = 1, dx = xgrid)$dC
  return(list(c(du1, du2)))
}


###################################################
### code chunk number 17: PDE.Rnw:387-389
###################################################
out <- ode.1D(func = wave, y = yini, times = times, parms = NULL,
         nspec = 2, method = "ode45", dimens = N, names = c("u", "v"))


###################################################
### code chunk number 18: wave
###################################################
matplot.1D(out, which = "u", subset = time %in% seq(0, 50, by = 10), 
   type = "l", col = c("black", rep("darkgrey", 5)), lwd = 2,
   grid = x, xlim = c(-50,50))
legend("topright", lty = 1:6, lwd = 2, col = c("black", rep("darkgrey", 5)), 
       legend = paste("t = ",seq(0, 50, by = 10)))


###################################################
### code chunk number 19: wave
###################################################
matplot.1D(out, which = "u", subset = time %in% seq(0, 50, by = 10), 
   type = "l", col = c("black", rep("darkgrey", 5)), lwd = 2,
   grid = x, xlim = c(-50,50))
legend("topright", lty = 1:6, lwd = 2, col = c("black", rep("darkgrey", 5)), 
       legend = paste("t = ",seq(0, 50, by = 10)))


###################################################
### code chunk number 20: waveimg
###################################################
par(mar=c(0,0,0,0))
image(out, which = "u", method = "persp", main = "",
        border = NA, col = "lightblue", box = FALSE, 
        shade = 0.5, theta = 0, phi = 60)                                                     


###################################################
### code chunk number 21: PDE.Rnw:465-473
###################################################
require(ReacTran)
pde2D <- function (t, y, parms) {
  CONC <- matrix(nr = n, nc = n, y)
  Tran <- tran.2D(CONC, D.x = Dx, D.y = Dy, dx = dx,  dy = dy)
  dCONC <- Tran$dC - r * CONC
  dCONC[ii]<- dCONC[ii] + p
  return(list(as.vector(dCONC)))
}


###################################################
### code chunk number 22: PDE.Rnw:482-488
###################################################
n  <- 100
dy <- dx <- 1
Dy <- Dx <- 2
r  <- 0.001
p  <- runif(50)
ii <- trunc(cbind(runif(50)*n, runif(50)*n) + 1)


###################################################
### code chunk number 23: PDE.Rnw:496-502
###################################################
require(rootSolve)
Conc0 <- matrix(nr = n, nc = n, 10.)
print(system.time(
  ST <- steady.2D(y = Conc0, func = pde2D, parms = NULL, dimens = c(n, n),
   lrw = 600000)
))


###################################################
### code chunk number 24: D2
###################################################
image(ST, main = "steady-state 2-D PDE")


###################################################
### code chunk number 25: D2fig
###################################################
image(ST, main = "steady-state 2-D PDE")


