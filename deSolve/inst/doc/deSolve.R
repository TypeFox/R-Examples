### R code from vignette source 'deSolve.Rnw'

###################################################
### code chunk number 1: preliminaries
###################################################
library("deSolve")
options(prompt = "> ")
options(width=70)


###################################################
### code chunk number 2: deSolve.Rnw:181-184
###################################################
parameters <- c(a = -8/3,
                b = -10,
                c =  28)


###################################################
### code chunk number 3: deSolve.Rnw:192-195
###################################################
state <- c(X = 1,
           Y = 1,
           Z = 1)


###################################################
### code chunk number 4: deSolve.Rnw:222-233
###################################################
Lorenz<-function(t, state, parameters) {
  with(as.list(c(state, parameters)),{
    # rate of change
    dX <- a*X + Y*Z
    dY <- b * (Y-Z)
    dZ <- -X*Y + c*Y - Z

    # return the rate of change
    list(c(dX, dY, dZ))
  })   # end with(as.list ...
}


###################################################
### code chunk number 5: deSolve.Rnw:243-244
###################################################
times <- seq(0, 100, by = 0.01)


###################################################
### code chunk number 6: deSolve.Rnw:259-262
###################################################
library(deSolve)
out <- ode(y = state, times = times, func = Lorenz, parms = parameters)
head(out)


###################################################
### code chunk number 7: ode
###################################################
par(oma = c(0, 0, 3, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "X"], out[, "Z"], pch = ".")

mtext(outer = TRUE, side = 3, "Lorenz model", cex = 1.5)


###################################################
### code chunk number 8: figode
###################################################
par(oma = c(0, 0, 3, 0))
plot(out, xlab = "time", ylab = "-")
plot(out[, "X"], out[, "Z"], pch = ".")

mtext(outer = TRUE, side = 3, "Lorenz model", cex = 1.5)


###################################################
### code chunk number 9: deSolve.Rnw:316-319
###################################################
outb <- radau(state, times, Lorenz, parameters, atol = 1e-4, rtol = 1e-4)
outc <- ode(state, times, Lorenz, parameters, method = "radau",
             atol = 1e-4, rtol = 1e-4)


###################################################
### code chunk number 10: deSolve.Rnw:335-341
###################################################
print(system.time(out1 <- rk4   (state, times, Lorenz, parameters)))
print(system.time(out2 <- lsode (state, times, Lorenz, parameters)))
print(system.time(out  <- lsoda (state, times, Lorenz, parameters)))
print(system.time(out  <- lsodes(state, times, Lorenz, parameters)))
print(system.time(out  <- daspk (state, times, Lorenz, parameters)))
print(system.time(out  <- vode  (state, times, Lorenz, parameters)))


###################################################
### code chunk number 11: deSolve.Rnw:359-360
###################################################
rkMethod()


###################################################
### code chunk number 12: deSolve.Rnw:369-370
###################################################
rkMethod("rk23")


###################################################
### code chunk number 13: deSolve.Rnw:383-404
###################################################
func <- function(t, x, parms) {
  with(as.list(c(parms, x)),{
    dP  <- a * P      - b * C * P
    dC  <- b * P * C  - c * C
    res <- c(dP, dC)
    list(res)
  })
}

rKnew <- rkMethod(ID = "midpoint",
  varstep = FALSE,
  A       = c(0, 1/2),
  b1      = c(0, 1),
  c       = c(0, 1/2),
  stage   = 2,
  Qerr    = 1
)

out <- ode(y = c(P = 2, C = 1), times = 0:100, func,
      parms = c(a = 0.1, b = 0.1, c = 0.1), method = rKnew)
head(out)


###################################################
### code chunk number 14: deSolve.Rnw:438-440
###################################################
diagnostics(out1)
diagnostics(out2)


###################################################
### code chunk number 15: deSolve.Rnw:444-445
###################################################
summary(out1)


###################################################
### code chunk number 16: deSolve.Rnw:519-527
###################################################
Aphid <- function(t, APHIDS, parameters) {
  deltax     <- c (0.5, rep(1, numboxes - 1), 0.5)
  Flux       <- -D * diff(c(0, APHIDS, 0)) / deltax
  dAPHIDS    <- -diff(Flux) / delx + APHIDS * r

  # the return value
  list(dAPHIDS )
} # end


###################################################
### code chunk number 17: deSolve.Rnw:532-539
###################################################
D         <- 0.3    # m2/day  diffusion rate
r         <- 0.01   # /day    net growth rate
delx      <- 1      # m       thickness of boxes
numboxes  <- 60

# distance of boxes on plant, m, 1 m intervals
Distance  <- seq(from = 0.5, by = delx, length.out = numboxes)


###################################################
### code chunk number 18: deSolve.Rnw:544-548
###################################################
# Initial conditions:  # ind/m2
APHIDS        <- rep(0, times = numboxes)
APHIDS[30:31] <- 1
state         <- c(APHIDS = APHIDS)      # initialise state variables


###################################################
### code chunk number 19: deSolve.Rnw:555-559
###################################################
times <-seq(0, 200, by = 1)
print(system.time(
  out <- ode.1D(state, times, Aphid, parms = 0, nspec = 1, names = "Aphid")
))


###################################################
### code chunk number 20: deSolve.Rnw:565-566
###################################################
head(out[,1:5])


###################################################
### code chunk number 21: deSolve.Rnw:570-571
###################################################
summary(out)


###################################################
### code chunk number 22: deSolve.Rnw:604-606
###################################################
data <- cbind(dist  = c(0,10, 20,  30,  40, 50, 60),
              Aphid = c(0,0.1,0.25,0.5,0.25,0.1,0))


###################################################
### code chunk number 23: matplot1d
###################################################
par (mfrow = c(1,2))
matplot.1D(out, grid = Distance, type = "l", mfrow = NULL,
    subset = time %in% seq(0, 200, by = 10),
   obs = data, obspar = list(pch = 18, cex = 2, col="red"))
plot.1D(out, grid = Distance, type = "l", mfrow = NULL,
    subset = time == 100,
   obs = data, obspar = list(pch = 18, cex = 2, col="red"))


###################################################
### code chunk number 24: matplot1d
###################################################
par (mfrow = c(1,2))
matplot.1D(out, grid = Distance, type = "l", mfrow = NULL,
    subset = time %in% seq(0, 200, by = 10),
   obs = data, obspar = list(pch = 18, cex = 2, col="red"))
plot.1D(out, grid = Distance, type = "l", mfrow = NULL,
    subset = time == 100,
   obs = data, obspar = list(pch = 18, cex = 2, col="red"))


###################################################
### code chunk number 25: deSolve.Rnw:672-687
###################################################
daefun <- function(t, y, dy, parameters) {
  res1 <- dy[1] + y[1] - y[2]
  res2 <- y[2] * y[1] - t

  list(c(res1, res2))
}

library(deSolve)
yini  <- c(1, 0)
dyini <- c(1, 0)
times <- seq(0, 10, 0.1)

## solver
system.time(out <- daspk(y = yini, dy = dyini,
                         times = times, res = daefun, parms = 0))


###################################################
### code chunk number 26: dae
###################################################
matplot(out[,1], out[,2:3], type = "l", lwd = 2,
        main = "dae", xlab = "time", ylab = "y")


###################################################
### code chunk number 27: figdae
###################################################
matplot(out[,1], out[,2:3], type = "l", lwd = 2,
        main = "dae", xlab = "time", ylab = "y")


###################################################
### code chunk number 28: deSolve.Rnw:720-730
###################################################
pendulum <- function (t, Y, parms) {
  with (as.list(Y),
    list(c(u,
           v,
          -lam * x,
          -lam * y - 9.8,
           x^2 + y^2 -1
     ))
  )
}


###################################################
### code chunk number 29: deSolve.Rnw:733-734
###################################################
yini <- c(x = 1, y = 0, u = 0, v = 1, lam = 1)


###################################################
### code chunk number 30: deSolve.Rnw:737-740
###################################################
M <- diag(nrow = 5)
M[5, 5] <- 0
M


###################################################
### code chunk number 31: deSolve.Rnw:744-748
###################################################
index <- c(2, 2, 1)
times <- seq(from = 0, to = 10, by = 0.01)
out  <- radau (y = yini, func = pendulum, parms = NULL,
               times = times, mass = M, nind = index)


###################################################
### code chunk number 32: pendulum
###################################################
plot(out, type = "l", lwd = 2)
plot(out[, c("x", "y")], type = "l", lwd = 2)


###################################################
### code chunk number 33: pendulum
###################################################
plot(out, type = "l", lwd = 2)
plot(out[, c("x", "y")], type = "l", lwd = 2)


###################################################
### code chunk number 34: deSolve.Rnw:782-795
###################################################
ZODE2 <- function(Time, State, Pars) {
  with(as.list(State), {
       df <- 1i * f
       dg <- -1i * g * g * f
       return(list(c(df, dg)))
  })
}

yini  <- c(f = 1+0i, g = 1/2.1+0i)
times <- seq(0, 2 * pi, length = 100)

out   <- zvode(func = ZODE2, y = yini, parms = NULL, times = times,
  atol = 1e-10, rtol = 1e-10)


###################################################
### code chunk number 35: deSolve.Rnw:807-809
###################################################
analytical <- cbind(f = exp(1i*times), g = 1/(exp(1i*times)+1.1))
tail(cbind(out[,2], analytical[,1]))


###################################################
### code chunk number 36: deSolve.Rnw:822-833
###################################################
f1 <- function  (t, y, parms) {
  ydot <- vector(len = 5)

  ydot[1] <-  0.1*y[1] -0.2*y[2]
  ydot[2] <- -0.3*y[1] +0.1*y[2] -0.2*y[3]
  ydot[3] <-           -0.3*y[2] +0.1*y[3] -0.2*y[4]
  ydot[4] <-                     -0.3*y[3] +0.1*y[4] -0.2*y[5]
  ydot[5] <-                               -0.3*y[4] +0.1*y[5]

  return(list(ydot))
}


###################################################
### code chunk number 37: deSolve.Rnw:838-840
###################################################
yini  <- 1:5
times <- 1:20


###################################################
### code chunk number 38: deSolve.Rnw:847-848
###################################################
out   <- lsode(yini, times, f1, parms = 0, jactype = "fullint")


###################################################
### code chunk number 39: deSolve.Rnw:855-864
###################################################
fulljac <- function  (t, y, parms) {
  jac <- matrix(nrow = 5, ncol = 5, byrow = TRUE,
                data = c(0.1, -0.2,  0  ,  0  ,  0  ,
                        -0.3,  0.1, -0.2,  0  ,  0  ,
                         0  , -0.3,  0.1, -0.2,  0  ,
                         0  ,  0  , -0.3,  0.1, -0.2,
                         0  ,  0  ,  0  , -0.3,  0.1))
  return(jac)
}


###################################################
### code chunk number 40: deSolve.Rnw:869-871
###################################################
out2 <- lsode(yini, times, f1, parms = 0, jactype = "fullusr",
              jacfunc = fulljac)


###################################################
### code chunk number 41: deSolve.Rnw:878-880
###################################################
out3 <- lsode(yini, times, f1, parms = 0, jactype = "bandint",
              bandup = 1, banddown = 1)


###################################################
### code chunk number 42: deSolve.Rnw:885-892
###################################################
bandjac <- function  (t, y, parms) {
  jac <- matrix(nrow = 3, ncol = 5, byrow = TRUE,
                data = c( 0  , -0.2, -0.2, -0.2, -0.2,
                          0.1,  0.1,  0.1,  0.1,  0.1,
                         -0.3, -0.3, -0.3, -0.3,  0))
  return(jac)
}


###################################################
### code chunk number 43: deSolve.Rnw:897-899
###################################################
out4 <- lsode(yini, times, f1, parms = 0, jactype = "bandusr",
              jacfunc = bandjac, bandup = 1, banddown = 1)


###################################################
### code chunk number 44: deSolve.Rnw:905-906
###################################################
out5  <- lsode(yini, times, f1, parms = 0, mf = 10)


###################################################
### code chunk number 45: deSolve.Rnw:937-943
###################################################
eventmod <- function(t, var, parms) {
  list(dvar = -0.1*var)
}

yini <- c(v1 = 1, v2 = 2)
times <- seq(0, 10, by = 0.1)


###################################################
### code chunk number 46: deSolve.Rnw:950-954
###################################################
eventdat <- data.frame(var = c("v1", "v2", "v2", "v1"), time = c(1, 1, 5, 9),
  value = c(1, 2, 3, 4), method = c("add", "mult", "rep", "add"))

eventdat


###################################################
### code chunk number 47: deSolve.Rnw:959-961
###################################################
out <- ode(func = eventmod, y = yini, times = times, parms = NULL,
  events = list(data = eventdat))


###################################################
### code chunk number 48: event1
###################################################
plot(out, type = "l", lwd = 2)


###################################################
### code chunk number 49: figevent1
###################################################
plot(out, type = "l", lwd = 2)


###################################################
### code chunk number 50: deSolve.Rnw:983-988
###################################################
ballode<- function(t, y, parms) {
  dy1 <- y[2]
  dy2 <- -9.8
  list(c(dy1, dy2))
}


###################################################
### code chunk number 51: deSolve.Rnw:995-996
###################################################
root <- function(t, y, parms) y[1]


###################################################
### code chunk number 52: deSolve.Rnw:1001-1006
###################################################
event <- function(t, y, parms) {
 y[1]<- 0
 y[2]<- -0.9 * y[2]
 return(y)
}


###################################################
### code chunk number 53: deSolve.Rnw:1012-1017
###################################################
yini  <- c(height = 0, v = 20)
times <- seq(from = 0, to = 20, by = 0.01)

out   <- lsode(times = times, y = yini, func = ballode, parms = NULL,
  events = list(func = event, root = TRUE), rootfun = root)


###################################################
### code chunk number 54: event2
###################################################
plot(out, which = "height", type = "l",lwd = 2,
  main = "bouncing ball", ylab = "height")


###################################################
### code chunk number 55: figevent2
###################################################
plot(out, which = "height", type = "l",lwd = 2,
  main = "bouncing ball", ylab = "height")


###################################################
### code chunk number 56: deSolve.Rnw:1066-1068
###################################################
times      <- seq(0, 1, 0.1)
eventtimes <- c(0.7, 0.9)


###################################################
### code chunk number 57: deSolve.Rnw:1073-1074
###################################################
eventtimes %in% times


###################################################
### code chunk number 58: deSolve.Rnw:1081-1083
###################################################
times2 <- round(times, 1)
times - times2


###################################################
### code chunk number 59: deSolve.Rnw:1094-1095
###################################################
eventtimes %in% times2


###################################################
### code chunk number 60: deSolve.Rnw:1100-1101
###################################################
all(eventtimes %in% times2)


###################################################
### code chunk number 61: deSolve.Rnw:1111-1114
###################################################
times <- 1:10
eventtimes <- c(1.3, 3.4, 4, 7.9, 8.5)
newtimes <- sort(unique(c(times, eventtimes)))


###################################################
### code chunk number 62: deSolve.Rnw:1120-1123
###################################################
times <- 1:10
eventtimes <- c(1.3, 3.4, 4, 7.9999999999999999, 8.5)
newtimes <- sort(c(eventtimes, cleanEventTimes(times, eventtimes)))


###################################################
### code chunk number 63: deSolve.Rnw:1152-1186
###################################################
library(deSolve)

#-----------------------------
# the derivative function
#-----------------------------
derivs <- function(t, y, parms) {
  if (t < 0)
    lag <- 19
  else
    lag <- lagvalue(t - 0.74)

  dy <- r * y * (1 - lag/m)
  list(dy, dy = dy)
}

#-----------------------------
# parameters
#-----------------------------

r <- 3.5; m <- 19

#-----------------------------
# initial values and times
#-----------------------------

yinit <- c(y = 19.001)
times <- seq(-0.74, 40, by = 0.01)

#-----------------------------
# solve the model
#-----------------------------

yout <- dede(y = yinit, times = times, func = derivs,
             parms = NULL, atol = 1e-10)


###################################################
### code chunk number 64: dde
###################################################
plot(yout, which = 1, type = "l", lwd = 2,
     main = "Lemming model", mfrow = c(1,2))
plot(yout[,2], yout[,3], xlab = "y", ylab = "dy", type = "l", lwd = 2)


###################################################
### code chunk number 65: figdde
###################################################
plot(yout, which = 1, type = "l", lwd = 2,
     main = "Lemming model", mfrow = c(1,2))
plot(yout[,2], yout[,3], xlab = "y", ylab = "dy", type = "l", lwd = 2)


###################################################
### code chunk number 66: deSolve.Rnw:1223-1235
###################################################
Stages <- c("DS 1yr", "DS 2yr", "R small", "R medium", "R large", "F")

NumStages   <- length(Stages)

# Population matrix
A <- matrix(nrow = NumStages, ncol = NumStages, byrow = TRUE, data = c(
                   0,      0,      0,      0,      0,      322.38,
                   0.966,  0,      0,      0,      0,      0     ,
                   0.013,  0.01,   0.125,  0,      0,      3.448 ,
                   0.007,  0,      0.125,  0.238,  0,      30.170,
                   0.008,  0,      0.038,  0.245,  0.167,  0.862 ,
                   0,      0,      0,      0.023,  0.75,   0      )  )


###################################################
### code chunk number 67: deSolve.Rnw:1240-1244
###################################################
Teasel <- function (t, y, p) {
  yNew <-  A %*% y
  list (yNew / sum(yNew))
}


###################################################
### code chunk number 68: deSolve.Rnw:1247-1249
###################################################
out <- ode(func = Teasel, y = c(1, rep(0, 5) ), times = 0:50,
           parms = 0, method = "iteration")


###################################################
### code chunk number 69: difference
###################################################
matplot(out[,1], out[,-1], main = "Teasel stage distribution", type = "l")
legend("topright", legend = Stages, lty = 1:6, col = 1:6)


###################################################
### code chunk number 70: difference
###################################################
matplot(out[,1], out[,-1], main = "Teasel stage distribution", type = "l")
legend("topright", legend = Stages, lty = 1:6, col = 1:6)


###################################################
### code chunk number 71: deSolve.Rnw:1292-1296
###################################################
library(deSolve)

combustion <- function (t, y, parms)
  list(y^2 * (1-y) )


###################################################
### code chunk number 72: deSolve.Rnw:1298-1300
###################################################
yini  <- 0.01
times <- 0 : 200


###################################################
### code chunk number 73: deSolve.Rnw:1302-1306
###################################################
out  <- ode(times = times, y = yini,   parms = 0, func = combustion)
out2 <- ode(times = times, y = yini*2, parms = 0, func = combustion)
out3 <- ode(times = times, y = yini*3, parms = 0, func = combustion)
out4 <- ode(times = times, y = yini*4, parms = 0, func = combustion)


###################################################
### code chunk number 74: plotdeSolve
###################################################
plot(out, out2, out3, out4, main = "combustion")
legend("bottomright", lty = 1:4, col = 1:4, legend = 1:4, title = "yini*i")


###################################################
### code chunk number 75: plotdeSolve
###################################################
plot(out, out2, out3, out4, main = "combustion")
legend("bottomright", lty = 1:4, col = 1:4, legend = 1:4, title = "yini*i")


###################################################
### code chunk number 76: deSolve.Rnw:1336-1337
###################################################
head(ccl4data)


###################################################
### code chunk number 77: deSolve.Rnw:1340-1343
###################################################
obs <- subset (ccl4data, animal == "A", c(time, ChamberConc))
names(obs) <- c("time", "CP")
head(obs)


###################################################
### code chunk number 78: deSolve.Rnw:1349-1363
###################################################
parms <- c(0.182, 4.0, 4.0, 0.08, 0.04, 0.74, 0.05, 0.15, 0.32, 16.17,
            281.48, 13.3, 16.17, 5.487, 153.8, 0.04321671,
            0.40272550, 951.46, 0.02, 1.0,  3.80000000)
yini <- c(AI = 21, AAM = 0, AT = 0, AF = 0, AL = 0, CLT = 0, AM = 0)

out <- ccl4model(times = seq(0, 6, by = 0.05), y = yini, parms = parms)

par2 <- parms
par2[1] <- 0.1
out2 <- ccl4model(times = seq(0, 6, by = 0.05), y = yini, parms = par2)

par3 <- parms
par3[1] <- 0.05
out3 <- ccl4model(times = seq(0, 6, by = 0.05), y = yini, parms = par3)


###################################################
### code chunk number 79: plotobs
###################################################
plot(out, out2, out3, which = c("AI", "MASS", "CP"),
     col = c("black", "red", "green"), lwd = 2,
     obs = obs, obspar = list(pch = 18, col = "blue", cex = 1.2))
legend("topright", lty = c(1,2,3,NA), pch = c(NA, NA, NA, 18),
  col = c("black", "red", "green", "blue"), lwd = 2,
  legend = c("par1", "par2", "par3", "obs"))


###################################################
### code chunk number 80: plotobs
###################################################
plot(out, out2, out3, which = c("AI", "MASS", "CP"),
     col = c("black", "red", "green"), lwd = 2,
     obs = obs, obspar = list(pch = 18, col = "blue", cex = 1.2))
legend("topright", lty = c(1,2,3,NA), pch = c(NA, NA, NA, 18),
  col = c("black", "red", "green", "blue"), lwd = 2,
  legend = c("par1", "par2", "par3", "obs"))


###################################################
### code chunk number 81: deSolve.Rnw:1389-1391
###################################################
obs2 <- data.frame(time = 6, MASS = 12)
obs2


###################################################
### code chunk number 82: obs2
###################################################
plot(out, out2, out3, lwd = 2,
     obs = list(obs, obs2),
     obspar = list(pch = c(16, 18), col = c("blue", "black"),
                   cex = c(1.2 , 2))
     )


###################################################
### code chunk number 83: plotobs2
###################################################
plot(out, out2, out3, lwd = 2,
     obs = list(obs, obs2),
     obspar = list(pch = c(16, 18), col = c("blue", "black"),
                   cex = c(1.2 , 2))
     )


###################################################
### code chunk number 84: hist
###################################################
hist(out, col = grey(seq(0, 1, by = 0.1)), mfrow = c(3, 4))


###################################################
### code chunk number 85: plothist
###################################################
hist(out, col = grey(seq(0, 1, by = 0.1)), mfrow = c(3, 4))


###################################################
### code chunk number 86: deSolve.Rnw:1450-1452
###################################################
options(prompt = " ")
options(continue = " ")


###################################################
### code chunk number 87: deSolve.Rnw:1455-1479
###################################################
lvmod <- function (time, state, parms, N, rr, ri, dr, dri) {
  with (as.list(parms), {
    PREY <- state[1:N]
    PRED <- state[(N+1):(2*N)]

    ## Fluxes due to diffusion
    ## at internal and external boundaries: zero gradient
    FluxPrey <- -Da * diff(c(PREY[1], PREY, PREY[N]))/dri
    FluxPred <- -Da * diff(c(PRED[1], PRED, PRED[N]))/dri

    ## Biology: Lotka-Volterra model
    Ingestion     <- rIng  * PREY * PRED
    GrowthPrey    <- rGrow * PREY * (1-PREY/cap)
    MortPredator  <- rMort * PRED

    ## Rate of change = Flux gradient + Biology
    dPREY    <- -diff(ri * FluxPrey)/rr/dr   +
                GrowthPrey - Ingestion
    dPRED    <- -diff(ri * FluxPred)/rr/dr   +
                Ingestion * assEff - MortPredator

    return (list(c(dPREY, dPRED)))
  })
}


###################################################
### code chunk number 88: deSolve.Rnw:1481-1483
###################################################
options(prompt = " ")
options(continue = " ")


###################################################
### code chunk number 89: deSolve.Rnw:1486-1500
###################################################
R  <- 20                        # total radius of surface, m
N  <- 100                       # 100 concentric circles
dr <- R/N                       # thickness of each layer
r  <- seq(dr/2,by = dr,len = N) # distance of center to mid-layer
ri <- seq(0,by = dr,len = N+1)  # distance to layer interface
dri <- dr                       # dispersion distances

parms <- c(Da     = 0.05,       # m2/d, dispersion coefficient
           rIng   = 0.2,        # /day, rate of ingestion
           rGrow  = 1.0,        # /day, growth rate of prey
           rMort  = 0.2 ,       # /day, mortality rate of pred
           assEff = 0.5,        # -, assimilation efficiency
           cap    = 10)         # density, carrying capacity



###################################################
### code chunk number 90: deSolve.Rnw:1503-1513
###################################################
state    <- rep(0, 2 * N)
state[1] <- state[N + 1] <- 10

times  <- seq(0, 200, by = 1)   # output wanted at these time intervals

print(system.time(
  out <- ode.1D(y = state, times = times, func = lvmod, parms = parms,
                nspec = 2, names = c("PREY", "PRED"),
                N = N, rr = r, ri = ri, dr = dr, dri = dri)
))


###################################################
### code chunk number 91: deSolve.Rnw:1516-1517
###################################################
summary(out)


###################################################
### code chunk number 92: deSolve.Rnw:1521-1523
###################################################
p10 <- subset(out, select = "PREY", subset = time == 10)
head(p10, n = 5)


###################################################
### code chunk number 93: deSolve.Rnw:1569-1574
###################################################
Simple2D <- function(t, Y, par) {
  y  <- matrix(nrow = nx, ncol = ny, data = Y)  # vector to 2-D matrix
  dY <- - r_x2y2 * y                            # consumption
  return(list(dY))
}


###################################################
### code chunk number 94: deSolve.Rnw:1578-1585
###################################################
dy <- dx <- 1  # grid size
nx <- ny <- 100
x  <- seq (dx/2, by = dx, len = nx)
y  <- seq (dy/2, by = dy, len = ny)

# in each grid cell: consumption depending on position
r_x2y2 <- outer(x, y, FUN=function(x,y) ((x-50)^2 + (y-50)^2)*1e-4)


###################################################
### code chunk number 95: deSolve.Rnw:1589-1592
###################################################
C <- matrix(nrow = nx, ncol = ny, 1)
ODE3 <- ode.2D(y = C, times = 1:100, func = Simple2D, parms = NULL,
               dimens = c(nx, ny), names = "C", method = "ode45")


###################################################
### code chunk number 96: deSolve.Rnw:1595-1598
###################################################
summary(ODE3)
t50 <-  matrix(nrow = nx, ncol = ny,
     data = subset(ODE3, select = "C", subset = (time == 50)))


###################################################
### code chunk number 97: twoD
###################################################
par(mfrow = c(1, 2))
contour(x, y,  r_x2y2, main = "consumption")
contour(x, y, t50, main = "Y(t = 50)")


###################################################
### code chunk number 98: twoD
###################################################
par(mfrow = c(1, 2))
contour(x, y,  r_x2y2, main = "consumption")
contour(x, y, t50, main = "Y(t = 50)")


###################################################
### code chunk number 99: deSolve.Rnw:1628-1636
###################################################
PCmod <- function(t, x, parms)  {
  with(as.list(c(parms, x)), {
    dP <- c*P   - d*C*P      # producer
    dC <- e*P*C - f*C        # consumer
    res <- c(dP, dC)
    list(res)
  })
}


###################################################
### code chunk number 100: deSolve.Rnw:1643-1644
###################################################
parms  <- c(c = 10, d = 0.1, e = 0.1, f = 0.1)


###################################################
### code chunk number 101: deSolve.Rnw:1650-1656
###################################################
xstart <- c(P = 0.5, C = 1)
times  <- seq(0, 200, 0.1)

out <- ode(y = xstart, times = times,
  func = PCmod, parms = parms)
tail(out)


###################################################
### code chunk number 102: deSolve.Rnw:1677-1681
###################################################
out <- ode(y = xstart,times = times, func = PCmod,
                         parms = parms, atol = 0)
matplot(out[,1], out[,2:3], type = "l",
        xlab = "time", ylab = "Producer, Consumer")


###################################################
### code chunk number 103: deSolve.Rnw:1737-1761
###################################################
LVmod <- function(Time, State, Pars) {
  with(as.list(c(State, Pars)), {
    Ingestion    <- rIng  * Prey * Predator
    GrowthPrey   <- rGrow * Prey * (1 - Prey/K)
    MortPredator <- rMort * Predator

    dPrey        <- GrowthPrey - Ingestion
    dPredator    <- Ingestion * assEff - MortPredator

    return(list(c(dPrey, dPredator)))
  })
}

pars    <- c(rIng   = 0.2,    # /day, rate of ingestion
             rGrow  = 1.0,    # /day, growth rate of prey
             rMort  = 0.2 ,   # /day, mortality rate of predator
             assEff = 0.5,    # -, assimilation efficiency
             K      = 10)     # mmol/m3, carrying capacity

yini    <- c(Prey = 1, Predator = 2)
times   <- seq(0, 200, by = 1)
out     <- ode(func = LVmod, y = yini,
              parms = pars, times = times)



###################################################
### code chunk number 104: deSolve.Rnw:1773-1776
###################################################
pars["rIng"] <- 100
out2 <- ode(func = LVmod, y = yini,
           parms = pars, times = times)


###################################################
### code chunk number 105: err
###################################################
plot(out2, type = "l", lwd = 2, main = "corrupt Lotka-Volterra model")


###################################################
### code chunk number 106: deSolve.Rnw:1825-1828
###################################################
pars["rIng"] <- 100
out3 <- ode(func = LVmod, y = yini, parms = pars,
          times = times, method = "ode45", atol = 1e-14, rtol = 1e-14)


###################################################
### code chunk number 107: err
###################################################
plot(out2, type = "l", lwd = 2, main = "corrupt Lotka-Volterra model")


