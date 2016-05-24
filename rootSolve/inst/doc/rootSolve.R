### R code from vignette source 'rootSolve.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("rootSolve")
options(prompt = "> ")
options(width=70)



###################################################
### code chunk number 2: fig1plot
###################################################
  fun <- function (x) cos(2*x)^3
  curve(fun(x), 0, 8)
  abline(h = 0, lty = 3)
  uni <- uniroot(fun, c(0, 8))$root
  points(uni, 0, pch = 16, cex = 2)


###################################################
### code chunk number 3: fig1
###################################################
  fun <- function (x) cos(2*x)^3
  curve(fun(x), 0, 8)
  abline(h = 0, lty = 3)
  uni <- uniroot(fun, c(0, 8))$root
  points(uni, 0, pch = 16, cex = 2)


###################################################
### code chunk number 4: uniall
###################################################
  curve(fun(x), 0, 8)
  abline(h = 0, lty = 3)
  All <- uniroot.all(fun, c(0, 8))
  points(All, y = rep(0, length(All)), pch = 16, cex = 2)


###################################################
### code chunk number 5: fig2
###################################################
  curve(fun(x), 0, 8)
  abline(h = 0, lty = 3)
  All <- uniroot.all(fun, c(0, 8))
  points(All, y = rep(0, length(All)), pch = 16, cex = 2)


###################################################
### code chunk number 6: rootSolve.Rnw:202-215
###################################################
model <- function(x) {
  F1 <- x[1]   + x[2]   + x[3]^2 -12
  F2 <- x[1]^2 - x[2]   + x[3] -2
  F3 <- 2*x[1] - x[2]^2 + x[3] -1
c(F1 = F1, F2 = F2, F3 = F3)
}

# first solution
(ss <- multiroot(f = model, start = c(1, 1, 1)))

# second solution; use different start values
(ss2 <- multiroot(model, c(0, 0, 0)))$root
model(ss2$root)   # the function value at the root


###################################################
### code chunk number 7: rootSolve.Rnw:234-244
###################################################
f2<-function(x)  {
 X <- matrix(nr = 5, x)
 X %*% X %*% X - matrix(nrow = 5, data = 1:25, byrow = TRUE)
}
print(system.time(
  x<-multiroot(f2, start = 1:25)$root
))
(X<-matrix(nrow = 5, x))

X%*%X%*%X


###################################################
### code chunk number 8: rootSolve.Rnw:262-266
###################################################
A <- matrix(nrow = 500, ncol = 500, runif(500*500))
B <- runif(500)

print(system.time(X1 <- solve(A, B)))


###################################################
### code chunk number 9: rootSolve.Rnw:273-274
###################################################
jfun <- function (x) A


###################################################
### code chunk number 10: rootSolve.Rnw:279-280
###################################################
fun <- function(x) A %*%x - B


###################################################
### code chunk number 11: rootSolve.Rnw:285-289
###################################################
print(system.time(
X <- multiroot(start = 1:500, f = fun, 
               jactype = "fullusr", jacfunc = jfun)
))


###################################################
### code chunk number 12: rootSolve.Rnw:292-293
###################################################
sum( (X1 - X$y)^2)


###################################################
### code chunk number 13: rootSolve.Rnw:359-373
###################################################
 model <- function(t, y, pars) {
  with (as.list(c(y, pars)),{

  oxicmin   = r*OM*(O2/(O2+ks))
  anoxicmin = r*OM*(1-O2/(O2+ks))* SO4/(SO4+ks2)

  dOM  = Flux - oxicmin - anoxicmin
  dO2  = -oxicmin       -2*rox*HS*(O2/(O2+ks)) + D*(BO2-O2)
  dSO4 = -0.5*anoxicmin   +rox*HS*(O2/(O2+ks)) + D*(BSO4-SO4)
  dHS  =  0.5*anoxicmin   -rox*HS*(O2/(O2+ks)) + D*(BHS-HS)

  list(c(dOM, dO2, dSO4, dHS), SumS = SO4+HS)
})
}


###################################################
### code chunk number 14: rootSolve.Rnw:378-388
###################################################
pars <- c(D = 1, Flux = 100, r = 0.1, rox = 1,
          ks = 1, ks2 = 1, BO2 = 100, BSO4 = 10000, BHS = 0)

y <- c(OM = 1, O2 = 1, SO4 = 1, HS = 1)

print(system.time(
  RS <- runsteady(y = y, fun = model, 
                  parms = pars, times = c(0, 1e5))
))
RS


###################################################
### code chunk number 15: rootSolve.Rnw:401-402
###################################################
stode(y = y, fun = model, parms = pars, pos = TRUE)


###################################################
### code chunk number 16: rootSolve.Rnw:441-449
###################################################
derivs <- function(t, y, parms, x, dx, N, y1, y6)  {

   d2y <- (c(y[-1],y6) -2*y + c(y1,y[-N])) /dx/dx
   dy  <- (c(y[-1],y6) - c(y1,y[-N])) /2/dx

   res <- d2y + dy/x + (1-1/(4*x*x))*y-sqrt(x)*cos(x)
   return(list(res))
}


###################################################
### code chunk number 17: rootSolve.Rnw:454-461
###################################################
dx     <- 0.001
x      <- seq(1, 6, by = dx)
N      <- length(x)
print(system.time(
y  <- steady.band(y = rep(1, N), time = 0, func = derivs, x = x,
                  dx = dx, N = N, y1 = 1, y6 = -0.5, nspec = 1)$y
))


###################################################
### code chunk number 18: steady1D
###################################################
plot(x, y, type = "l", 
     main = "5001 nonlinear equations - banded Jacobian")

curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+
      0.740071*sin(x)/sqrt(x)+1/4*x^(3/2)*sin(x),add=TRUE,type="p")

legend("topright", pch = c(NA, 1), lty = c(1, NA),
       c("numeric", "analytic"))


###################################################
### code chunk number 19: steady1D
###################################################
plot(x, y, type = "l", 
     main = "5001 nonlinear equations - banded Jacobian")

curve(0.0588713*cos(x)/sqrt(x)+1/4*sqrt(x)*cos(x)+
      0.740071*sin(x)/sqrt(x)+1/4*x^(3/2)*sin(x),add=TRUE,type="p")

legend("topright", pch = c(NA, 1), lty = c(1, NA),
       c("numeric", "analytic"))


###################################################
### code chunk number 20: rootSolve.Rnw:505-520
###################################################
O2BOD <- function(t, state, pars) {
  BOD <- state[1:N]
  O2  <- state[(N+1):(2*N)]

  FluxBOD <-  v*c(BOD_0,BOD)  # fluxes due to water transport
  FluxO2  <-  v*c(O2_0,O2)

  BODrate <- r*BOD*O2/(O2+10)  # 1-st order consumption, Monod in oxygen

#rate of change = -flux gradient - consumption  + reaeration (O2)
  dBOD         <- -diff(FluxBOD)/dx  - BODrate
  dO2          <- -diff(FluxO2)/dx   - BODrate + p*(O2sat-O2)

  return(list(c(dBOD = dBOD, dO2 = dO2), BODrate = BODrate))
 }


###################################################
### code chunk number 21: rootSolve.Rnw:527-544
###################################################
dx      <- 10        # grid size, meters
v       <- 1e2       # velocity, m/day
r       <- 0.1       # /day, first-order decay of BOD
p       <- 0.1       # /day, air-sea exchange rate
O2sat   <- 300       # mmol/m3 saturated oxygen conc
O2_0    <- 50        # mmol/m3 riverine oxygen conc
BOD_0   <- 1500      # mmol/m3 riverine BOD concentration

x       <- seq(dx/2, 10000, by = dx)  # m, distance from river
N       <- length(x)

state <- c(rep(200, N), rep(200, N))    # initial guess of state variables:

print(system.time(
  out   <- steady.1D (y = state, func = O2BOD, parms = NULL, 
    nspec = 2, pos = TRUE)
))


###################################################
### code chunk number 22: BOD
###################################################
mf <- par(mfrow = c(2, 2))
plot(out, grid = x, xlab = "Distance from river", mfrow = NULL,
     ylab = "mmol/m3", main = c("Oxygen", "BOD"), type = "l")
plot(out, which = "BODrate", grid = x, mfrow = NULL,
     xlab = "Distance from river",
     ylab = "mmol/m3/d", main = "BOD decay rate", type = "l")
par(mfrow = mf)


###################################################
### code chunk number 23: figBOD
###################################################
mf <- par(mfrow = c(2, 2))
plot(out, grid = x, xlab = "Distance from river", mfrow = NULL,
     ylab = "mmol/m3", main = c("Oxygen", "BOD"), type = "l")
plot(out, which = "BODrate", grid = x, mfrow = NULL,
     xlab = "Distance from river",
     ylab = "mmol/m3/d", main = "BOD decay rate", type = "l")
par(mfrow = mf)


###################################################
### code chunk number 24: rootSolve.Rnw:614-635
###################################################
diffusion2D <- function(t, conc, par)  {
   Conc     <- matrix(nrow = n, ncol = n, data = conc)  # vector to 2-D matrix
   dConc    <- -r*Conc*Conc    # consumption
   BND   <- rep(1, n)           # boundary concentration

   # constant production in certain cells
   dConc[ii]<-   dConc[ii]+  p

   #diffusion in X-direction; boundaries=imposed concentration

   Flux  <- -Dx * rbind(rep(0, n), (Conc[2:n,]-Conc[1:(n-1),]),
                        rep(0, n) )/dx
   dConc <- dConc - (Flux[2:(n+1),] - Flux[1:n,])/dx

   #diffusion in Y-direction
   Flux  <- -Dy * cbind(rep(0, n), (Conc[,2:n]-Conc[,1:(n-1)]),
                        rep(0, n))/dy
   dConc <- dConc - (Flux[,2:(n+1)]-Flux[,1:n])/dy

   return(list(as.vector(dConc)))
  }


###################################################
### code chunk number 25: rootSolve.Rnw:639-647
###################################################
  # parameters
  dy    <- dx <- 1    # grid size
  Dy    <- Dx <- 1.5  # diffusion coeff, X- and Y-direction
  r     <- 0.01       # 2-nd-order consumption rate (/time)
  p     <- 20         # 0-th order production rate (CONC/t)
  n     <- 100
  # 10 random cells where substance is produced at rate p
  ii    <- trunc(cbind(runif(10)*n+1, runif(10)*n+1))


###################################################
### code chunk number 26: rootSolve.Rnw:656-663
###################################################
  Conc0 <- matrix(nrow = n, ncol = n, 10.)

  print(system.time(
  ST3  <- steady.2D(Conc0, func = diffusion2D, parms = NULL,
                    pos = TRUE, dimens = c(n, n), lrw = 1000000,
                    atol = 1e-10, rtol = 1e-10, ctol = 1e-10)
  ))


###################################################
### code chunk number 27: steady2D
###################################################
image(ST3, main = "2-D diffusion+production", xlab = "x", ylab = "y",
      legend = TRUE)


###################################################
### code chunk number 28: fig2D
###################################################
image(ST3, main = "2-D diffusion+production", xlab = "x", ylab = "y",
      legend = TRUE)


###################################################
### code chunk number 29: rootSolve.Rnw:685-739
###################################################
diffusion3D <- function(t, Y, par)   {
   yy    <- array(dim=c(n,n,n),data=Y)  # vector to 3-D array
   dY   <- -r*yy        # consumption
   BND   <- rep(1,n)   # boundary concentration
   for (i in 1:n) {
     y <- yy[i,,]

     #diffusion in X-direction; boundaries=imposed concentration
     Flux <- -Dy * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dy
     dY[i,,]   <- dY[i,,] - (Flux[2:(n+1),]-Flux[1:n,])/dy

     #diffusion in Y-direction
     Flux <- -Dz * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dz
     dY[i,,]    <- dY[i,,] - (Flux[,2:(n+1)]-Flux[,1:n])/dz
   }
   for (j in 1:n) {
     y <- yy[,j,]

     #diffusion in X-direction; boundaries=imposed concentration
     Flux <- -Dx * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dx
     dY[,j,]   <- dY[,j,] - (Flux[2:(n+1),]-Flux[1:n,])/dx

     #diffusion in Y-direction
     Flux <- -Dz * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dz
     dY[,j,]    <- dY[,j,] - (Flux[,2:(n+1)]-Flux[,1:n])/dz
   }
   for (k in 1:n) {
     y <- yy[,,k]

     #diffusion in X-direction; boundaries=imposed concentration
     Flux <- -Dx * rbind(y[1,]-BND,(y[2:n,]-y[1:(n-1),]),BND-y[n,])/dx
     dY[,,k]   <- dY[,,k] - (Flux[2:(n+1),]-Flux[1:n,])/dx

     #diffusion in Y-direction
     Flux <- -Dy * cbind(y[,1]-BND,(y[,2:n]-y[,1:(n-1)]),BND-y[,n])/dy
     dY[,,k]    <- dY[,,k] - (Flux[,2:(n+1)]-Flux[,1:n])/dy
   }
   return(list(as.vector(dY)))
}

  # parameters
  dy    <- dx <- dz <-1   # grid size
  Dy    <- Dx <- Dz <-1   # diffusion coeff, X- and Y-direction
  r     <- 0.025     # consumption rate

  n  <- 10
  y  <- array(dim=c(n, n, n), data = 10.)

  print(system.time(
  ST3 <- steady.3D(y, func =diffusion3D, parms = NULL, pos = TRUE,
                   dimens = c(n, n, n), lrw = 100000,
                   atol = 1e-10, rtol = 1e-10, ctol = 1e-10, 
                   verbose = TRUE)
  ))


###################################################
### code chunk number 30: steady3D
###################################################
  image(ST3, mfrow = c(2, 2), add.contour = TRUE, legend = TRUE,
        dimselect = list(x = c(1, 4, 8, 10)))


###################################################
### code chunk number 31: fig3D
###################################################
  image(ST3, mfrow = c(2, 2), add.contour = TRUE, legend = TRUE,
        dimselect = list(x = c(1, 4, 8, 10)))


###################################################
### code chunk number 32: rootSolve.Rnw:831-837
###################################################
bvp22 <- function (y, xi) {
  dy2 <- diff(diff(c(ya, y, yb))/dx)/dx
  dy  <- 0.5*(diff(c(ya, y))/dx + diff(c(y, yb))/dx)
  
  return(xi*dy2+dy+y^2)
}


###################################################
### code chunk number 33: rootSolve.Rnw:841-846
###################################################
dx <- 0.001
x  <- seq(0, 1, by = dx)
N  <- length(x)
ya <- 0
yb <- 0.5


###################################################
### code chunk number 34: rootSolve.Rnw:851-859
###################################################
print(system.time(
  Y1<- multiroot.1D(f = bvp22, start = runif(N), nspec = 1, xi = 0.1)
)*1000)
  Y2<- multiroot.1D(f = bvp22, start = runif(N), nspec = 1, xi = 0.05)
print(system.time(
  Y3<- multiroot.1D(f = bvp22, start = runif(N), nspec = 1, xi = 0.01)
)*1000)



###################################################
### code chunk number 35: bnd
###################################################
plot(x, Y3$root, type = "l", col = "green", lwd = 2,
     main = "bvp test problem 22" , ylab = "y")
lines(x, Y2$root, col = "red", lwd = 2)
lines(x, Y1$root, col = "blue", lwd = 2)


###################################################
### code chunk number 36: figbnd
###################################################
plot(x, Y3$root, type = "l", col = "green", lwd = 2,
     main = "bvp test problem 22" , ylab = "y")
lines(x, Y2$root, col = "red", lwd = 2)
lines(x, Y1$root, col = "blue", lwd = 2)


###################################################
### code chunk number 37: rootSolve.Rnw:1069-1079
###################################################
# the banana function
   fun <- function(x)  100*(x[2] - x[1]^2)^2 + (1 - x[1])^2
# the minimum
   mm  <-nlm(fun, p=c(0,0))$estimate
# the Hessian
   (Hes <- hessian(fun,mm))
# the gradient
   (grad <- gradient(fun,mm,centered=TRUE))
# Hessian and gradient can also be estimated by nlm:
   nlm(fun, p=c(0,0), hessian=TRUE)


###################################################
### code chunk number 38: rootSolve.Rnw:1082-1083
###################################################
   solve(Hes)


###################################################
### code chunk number 39: rootSolve.Rnw:1091-1101
###################################################
mod <- function (t=0,y, parms=NULL,...)
{
 dy1 <-   y[1] + 2*y[2]
 dy2 <- 3*y[1] + 4*y[2] + 5*y[3]
 dy3 <-          6*y[2] + 7*y[3] + 8*y[4]
 dy4 <-                   9*y[3] +10*y[4]
 return(as.list(c(dy1, dy2, dy3, dy4)))
}
jacobian.full(y = c(1, 2, 3, 4), func = mod)
jacobian.band(y = c(1, 2, 3, 4), func = mod)


