### R code from vignette source 'laGP.Rnw'

###################################################
### code chunk number 1: laGP.Rnw:87-94
###################################################
library("laGP")
library("MASS")
library("lhs")
library("akima")
library("tgp")
options(prompt="R> ", width=65)
set.seed(1)


###################################################
### code chunk number 2: laGP.Rnw:240-242
###################################################
X <- matrix(seq(0, 2 * pi,length = 6), ncol = 1)
Z <- sin(X)


###################################################
### code chunk number 3: laGP.Rnw:248-250
###################################################
gp <- newGP(X, Z, 2, 1e-6, dK = TRUE)
mleGP(gp, tmax=20)


###################################################
### code chunk number 4: laGP.Rnw:266-269
###################################################
XX <- matrix(seq(-1, 2 * pi + 1, length = 499), ncol = ncol(X))
p <- predGP(gp, XX)
deleteGP(gp)


###################################################
### code chunk number 5: laGP.Rnw:276-280
###################################################
library("mvtnorm")
N <- 100
ZZ <- rmvt(N, p$Sigma, p$df)
ZZ <- ZZ + t(matrix(rep(p$mean, N), ncol = N))


###################################################
### code chunk number 6: sin
###################################################
matplot(XX, t(ZZ), col = "gray", lwd = 0.5, lty = 1, type = "l",
	bty = "n", main = "simple sinusoidal example", xlab = "x", 
	ylab = "Y(x) | thetahat")
points(X, Z, pch = 19)


###################################################
### code chunk number 7: laGP.Rnw:485-488
###################################################
x <- seq(-2, 2, by = 0.02)
X <- as.matrix(expand.grid(x, x))
N <- nrow(X)


###################################################
### code chunk number 8: laGP.Rnw:502-510
###################################################
f2d <- function(x)
  {
    g <- function(z)
      return(exp( - (z - 1)^2) + exp( -0.8 * (z + 1)^2) 
        - 0.05 * sin(8 * (z + 0.1)))
    -g(x[,1]) * g(x[,2])
  }
Y <- f2d(X)


###################################################
### code chunk number 9: laGP.Rnw:515-518
###################################################
Xref <- matrix(c(-1.725, 1.725), nrow = 1)
p.mspe <- laGP(Xref, 6, 50, X, Y, d = 0.1, method="mspe")
p.alc <- laGP(Xref, 6, 50, X, Y, d = 0.1, method="alc")


###################################################
### code chunk number 10: lagp
###################################################
Xi <- rbind(X[p.mspe$Xi, ], X[p.alc$Xi, ])
plot(X[p.mspe$Xi, ], xlab = "x1", ylab = "x2", type = "n", 
  main = "comparing local designs", xlim = range(Xi[ ,1]), 
  ylim = range(Xi[ ,2]))
text(X[p.mspe$Xi, ], labels = 1:length(p.mspe$Xi), cex = 0.7)
text(X[p.alc$Xi, ], labels = 1:length(p.alc$Xi), cex = 0.7, col = 2)
points(Xref[1], Xref[2], pch=19, col=3)
legend("topright", c("mspe", "alc"), text.col = c(1, 2), bty="n")


###################################################
### code chunk number 11: laGP.Rnw:564-569
###################################################
p <- rbind(c(p.mspe$mean, p.mspe$s2, p.mspe$df),
  c(p.alc$mean, p.alc$s2, p.alc$df))
colnames(p) <- c("mean", "s2", "df")
rownames(p) <- c("mspe", "alc")
p


###################################################
### code chunk number 12: laGP.Rnw:573-575
###################################################
p.mspe$mle
p.alc$mle


###################################################
### code chunk number 13: laGP.Rnw:586-587
###################################################
c(p.mspe$time, p.alc$time)


###################################################
### code chunk number 14: laGP.Rnw:638-641
###################################################
xx <- seq(-1.97, 1.95, by = 0.04)
XX <- as.matrix(expand.grid(xx, xx))
YY <- f2d(XX)


###################################################
### code chunk number 15: laGP.Rnw:652-655
###################################################
nth <- as.numeric(Sys.getenv("OMP_NUM_THREADS"))
if(is.na(nth)) nth <- 2
print(nth)


###################################################
### code chunk number 16: laGP.Rnw:661-662
###################################################
P.alc <- aGP(X, Y, XX, omp.threads = nth, verb = 0)


###################################################
### code chunk number 17: aggp
###################################################
persp(xx, xx, -matrix(P.alc$mean, ncol = length(xx)), phi=45, theta=45,
      main = "", xlab = "x1", ylab = "x2", zlab = "yhat(x)")


###################################################
### code chunk number 18: aggp-slice
###################################################
med <- 0.51
zs <- XX[, 2] == med
sv <- sqrt(P.alc$var[zs])
r <- range(c(-P.alc$mean[zs] + 2 * sv, -P.alc$mean[zs] - 2 * sv))
plot(XX[zs,1], -P.alc$mean[zs], type="l", lwd = 2, ylim = r, xlab = "x1",
     ylab = "predicted & true response", bty = "n",
     main = "slice through surface")
lines(XX[zs, 1], -P.alc$mean[zs] + 2 * sv, col = 2, lty = 2, lwd = 2)
lines(XX[zs, 1], -P.alc$mean[zs] - 2 * sv, col = 2, lty = 2, lwd = 2)
lines(XX[zs, 1], YY[zs], col = 3, lwd = 2, lty = 3)


###################################################
### code chunk number 19: aggp-slice-ydiff
###################################################
diff <- P.alc$mean - YY
plot(XX[zs,1], diff[zs], type = "l", lwd = 2, 
     main = "systematic bias in prediction", 
     xlab = "x1", ylab = "y(x) - yhat(x)", bty = "n")


###################################################
### code chunk number 20: aggp-slice-d
###################################################
plot(XX[zs,1], P.alc$mle$d[zs], type = "l", lwd=2, 
     main = "spatially varying lengthscale", 
     xlab = "x1", ylab = "thetahat(x)", bty = "n")
df <- data.frame(y = log(P.alc$mle$d), XX)
lo <- loess(y ~ ., data = df, span = 0.01)
lines(XX[zs,1], exp(lo$fitted)[zs], col=2, lty=2, lwd=2)
legend("topright", "loess smoothed", col=2, lty=2, lwd=2, bty="n")


###################################################
### code chunk number 21: laGP.Rnw:823-824
###################################################
P.alc2 <- aGP(X, Y, XX, d = exp(lo$fitted), omp.threads = nth, verb = 0)


###################################################
### code chunk number 22: laGP.Rnw:831-834
###################################################
rmse <- data.frame(alc = sqrt(mean((P.alc$mean - YY)^2)), 
  alc2 = sqrt(mean((P.alc2$mean - YY)^2)))
rmse


###################################################
### code chunk number 23: laGP.Rnw:904-905
###################################################
p.alcray <- laGP(Xref, 6, 50, X, Y, d = 0.1, method = "alcray")


###################################################
### code chunk number 24: lagp-ray
###################################################
plot(X[p.alc$Xi,], xlab = "x1", ylab = "x2", type = "n", 
  main="comparing local designs", xlim = range(Xi[ ,1]), 
  ylim = range(Xi[ ,2]))
text(X[p.alc$Xi,], labels = 1:length(p.alc$Xi), cex = 0.7, col = 2)
text(X[p.alcray$Xi,], labels=1:length(p.mspe$Xi), cex=0.7, col = 3)
points(Xref[1], Xref[2], pch = 19, col = 3)
legend("topright", c("alc", "alcray"), text.col = c(2,3), bty = "n")


###################################################
### code chunk number 25: laGP.Rnw:928-929
###################################################
p.alcray$time


###################################################
### code chunk number 26: laGP.Rnw:933-936
###################################################
p <- rbind(p, c(p.alcray$mean, p.alcray$s2, p.alcray$df))
rownames(p)[3] <- c("alcray")
p


###################################################
### code chunk number 27: laGP.Rnw:945-950
###################################################
P.alcray <- aGP(X, Y, XX, method = "alcray", omp.threads = nth, verb = 0)
dfray <- data.frame(y = log(P.alcray$mle$d), XX)
loray <- loess(y ~ ., data = dfray, span = 0.01)
P.alcray2 <- aGP(X, Y, XX, method = "alcray", d = exp(loray$fitted), 
  omp.threads = nth, verb = 0)


###################################################
### code chunk number 28: laGP.Rnw:956-957
###################################################
c(P.alcray$time, P.alcray2$time)


###################################################
### code chunk number 29: laGP.Rnw:960-964
###################################################
rmse <- cbind(rmse, 
  data.frame(alcray=sqrt(mean((P.alcray$mean - YY)^2)), 
    alcray2=sqrt(mean((P.alcray2$mean - YY)^2))))
rmse


###################################################
### code chunk number 30: laGP.Rnw:1005-1019
###################################################
borehole <- function(x){
  rw <- x[1] * (0.15 - 0.05) + 0.05
  r <-  x[2] * (50000 - 100) + 100
  Tu <- x[3] * (115600 - 63070) + 63070
  Hu <- x[4] * (1110 - 990) + 990
  Tl <- x[5] * (116 - 63.1) + 63.1
  Hl <- x[6] * (820 - 700) + 700
  L <-  x[7] * (1680 - 1120) + 1120
  Kw <- x[8] * (12045 - 9855) + 9855
  m1 <- 2 * pi * Tu * (Hu - Hl)
  m2 <- log(r / rw)
  m3 <- 1 + 2 * L * Tu / (m2 * rw^2 * Kw) + Tu / Tl
  return(m1/m2/m3)
}


###################################################
### code chunk number 31: laGP.Rnw:1028-1032
###################################################
N <- 100000
Npred <- 1000
dim <- 8
library("lhs")


###################################################
### code chunk number 32: laGP.Rnw:1041-1047
###################################################
T <- 10
nas <- rep(NA, T)
times <- rmse <- data.frame(mspe = nas, mspe2 = nas, 
  alc.nomle = nas, alc = nas, alc2 = nas,
  nn.nomle = nas, nn=nas, big.nn.nomle = nas, big.nn = nas,
  big.alcray = nas, big.alcray2 = nas)


###################################################
### code chunk number 33: laGP.Rnw:1061-1124
###################################################
for(t in 1:T) {

  x <- randomLHS(N + Npred, dim)
  y <- apply(x, 1, borehole)
  ypred.0 <- y[-(1:N)]; y <- y[1:N]
  xpred <- x[-(1:N),]; x <- x[1:N,]

  formals(aGP)[c("omp.threads", "verb")] <- c(nth, 0)
  formals(aGP)[c("X", "Z", "XX")] <- list(x, y, xpred)

  out1<- aGP(d=list(mle = FALSE, start = 0.7))
  rmse$alc.nomle[t] <- sqrt(mean((out1$mean - ypred.0)^2))
  times$alc.nomle[t] <- out1$time
  
  out2 <- aGP(d = list(max = 20))
  rmse$alc[t] <- sqrt(mean((out2$mean - ypred.0)^2))
  times$alc[t] <- out2$time
  
  out3 <- aGP(d = list(start = out2$mle$d, max = 20))
  rmse$alc2[t] <- sqrt(mean((out3$mean - ypred.0)^2))
  times$alc2[t] <- out3$time

  out4 <- aGP(d = list(max = 20), method="alcray")
  rmse$alcray[t] <- sqrt(mean((out4$mean - ypred.0)^2))
  times$alcray[t] <- out4$time
  
  out5 <- aGP(d = list(start = out4$mle$d, max = 20), method="alcray")
  rmse$alcray2[t] <- sqrt(mean((out5$mean - ypred.0)^2))
  times$alcray2[t] <- out5$time

  out6<- aGP(d = list(max = 20), method="mspe")
  rmse$mspe[t] <- sqrt(mean((out6$mean - ypred.0)^2))
  times$mspe[t] <- out6$time
  
  out7 <- aGP(d = list(start = out6$mle$d, max = 20), method="mspe")
  rmse$mspe2[t] <- sqrt(mean((out7$mean - ypred.0)^2))
  times$mspe2[t] <- out7$time

  out8 <- aGP(d = list(mle = FALSE, start = 0.7), method = "nn")
  rmse$nn.nomle[t] <- sqrt(mean((out8$mean - ypred.0)^2))
  times$nn.nomle[t] <- out8$time

  out9 <- aGP(end = 200, d = list(mle = FALSE), method = "nn")
  rmse$big.nn.nomle[t] <- sqrt(mean((out9$mean - ypred.0)^2))
  times$big.nn.nomle[t] <- out9$time

  out10 <- aGP(d = list(max = 20), method = "nn")
  rmse$nn[t] <- sqrt(mean((out10$mean - ypred.0)^2))
  times$nn[t] <- out10$time

  out11 <- aGP(end = 200, d = list(max = 20), method="nn")
  rmse$big.nn[t] <- sqrt(mean((out11$mean - ypred.0)^2))
  times$big.nn[t] <- out11$time

  out12 <- aGP(end = 200, d = list(max = 20), method="alcray")
  rmse$big.alcray[t] <- sqrt(mean((out12$mean - ypred.0)^2))
  times$big.alcray[t] <- out12$time
  
  out13 <- aGP(end = 200, d = list(start = out12$mle$d, max = 20), 
    method="alcray")
  rmse$big.alcray2[t] <- sqrt(mean((out13$mean - ypred.0)^2))
  times$big.alcray2[t] <- out13 $time
}


###################################################
### code chunk number 34: laGP.Rnw:1132-1144
###################################################
timev <- apply(times, 2, mean, na.rm = TRUE)
rmsev <- apply(rmse, 2, mean)
tab <- cbind(timev, rmsev)
o <- order(rmsev, decreasing = FALSE)
tt <- rep(NA, length(rmsev))
for(i in 1:(length(o)-1)) {
  tto <- t.test(rmse[ ,o[i]], rmse[ ,o[i+1]], alternative = "less", 
    paired = TRUE)
  tt[o[i]] <- tto$p.value
}
tab <- cbind(tab, data.frame(tt))
tab[o, ]


###################################################
### code chunk number 35: laGP.Rnw:1209-1229
###################################################
thats <- matrix(NA, nrow = T, ncol = dim)
its <- rep(NA, T)
n <- 1000

g2 <- garg(list(mle = TRUE), y)
d2 <- darg(list(mle = TRUE, max = 100), x)

for(t in 1:T) {
  
  subs <- sample(1:N, n, replace = FALSE)

  gpsepi <- newGPsep(x[subs, ], y[subs], rep(d2$start, dim), g = 1/1000, 
    dK = TRUE)
  that <- mleGPsep(gpsepi, param = "d", tmin = d2$min, tmax = d2$max, 
    ab = d2$ab, maxit = 200)
  thats[t,] <- that$d
  its[t] <- that$its

  deleteGPsep(gpsepi)
}


###################################################
### code chunk number 36: thetas
###################################################
boxplot(thats, main = "distribution of thetas", xlab = "input", 
  ylab = "theta")


###################################################
### code chunk number 37: laGP.Rnw:1275-1281
###################################################
scales <- sqrt(apply(thats, 2, median))
xs <- x; xpreds <- xpred
for(j in 1:ncol(xs)) {
  xs[,j] <- xs[,j] / scales[j]
  xpreds[,j] <- xpreds[,j] / scales[j]
}


###################################################
### code chunk number 38: laGP.Rnw:1285-1286
###################################################
out14 <- aGP(xs, y, xpreds, d=list(start=1, max=20), method="alcray")


###################################################
### code chunk number 39: laGP.Rnw:1291-1292
###################################################
sqrt(mean((out14$mean - ypred.0)^2))


###################################################
### code chunk number 40: laGP.Rnw:1317-1324
###################################################
library("MASS")
d <- darg(NULL, mcycle[, 1, drop = FALSE])
g <- garg(list(mle = TRUE), mcycle[,2])
motogp <- newGP(mcycle[ , 1, drop=FALSE], mcycle[ ,2], d = d$start, 
  g = g$start, dK = TRUE)
jmleGP(motogp, drange = c(d$min, d$max), grange = c(d$min, d$max), 
  dab = d$ab, gab = g$ab)


###################################################
### code chunk number 41: laGP.Rnw:1329-1334
###################################################
XX <- matrix(seq(min(mcycle[ ,1]), max(mcycle[ ,1]), length = 100), 
  ncol = 1)
motogp.p <- predGP(motogp, XX = XX, lite = TRUE)
motoagp <- aGP(mcycle[ , 1, drop=FALSE], mcycle[,2], XX, end = 30, 
  d = d, g = g, verb = 0)


###################################################
### code chunk number 42: mcycle
###################################################
plot(mcycle, cex = 0.5, main = "motorcycle data")
lines(XX, motogp.p$mean, lwd = 2)
q1 <- qnorm(0.05, mean = motogp.p$mean, sd = sqrt(motogp.p$s2))
q2 <- qnorm(0.95, mean = motogp.p$mean, sd = sqrt(motogp.p$s2))
lines(XX, q1, lty = 2, lwd = 2)
lines(XX, q2, lty = 2, lwd = 2)
lines(XX, motoagp$mean, col = 2, lwd = 2)
q1 <- qnorm(0.05, mean = motoagp$mean, sd = sqrt(motoagp$var))
q2 <- qnorm(0.95, mean = motoagp$mean, sd = sqrt(motoagp$var))
lines(XX, q1, lty = 2, col = 2, lwd = 2)
lines(XX, q2, lty = 2, col = 2, lwd = 2)


###################################################
### code chunk number 43: laGP.Rnw:1376-1380
###################################################
X <- matrix(rep(mcycle[ ,1], 10), ncol = 1)
X <- X + rnorm(nrow(X), sd = 1)
Z <- rep(mcycle[ ,2], 10)
motoagp2 <- aGP(X, Z, XX, end = 30, d = d, g = g, verb = 0)


###################################################
### code chunk number 44: mcycle-rep
###################################################
plot(X, Z, main = "simulating a larger data setup", xlab = "times", 
  ylab = "accel")
lines(XX, motoagp2$mean, col = 2, lwd = 2)
q1 <- qnorm(0.05, mean = motoagp2$mean, sd = sqrt(motoagp2$var))
q2 <- qnorm(0.95, mean = motoagp2$mean, sd = sqrt(motoagp2$var))
lines(XX, q1, col = 2, lty = 2, lwd = 2)
lines(XX, q2, col = 2, lty = 2, lwd = 2)


###################################################
### code chunk number 45: laGP.Rnw:1628-1638
###################################################
M <- function(x,u) 
  {
    x <- as.matrix(x)
    u <- as.matrix(u)
    out <- (1 - exp(-1 / (2 * x[,2]))) 
    out <- out * (1000 * u[,1] * x[,1]^3 + 1900 * x[ ,1]^2 
      + 2092 * x[ ,1] + 60) 
    out <- out / (100 * u[,2] * x[,1]^3 + 500 * x[ ,1]^2 + 4 * x[ ,1] + 20)  
    return(out)
  }


###################################################
### code chunk number 46: laGP.Rnw:1642-1648
###################################################
bias <- function(x) 
  {
    x <- as.matrix(x)   
    out <- 2 * (10 * x[ ,1]^2 + 4 * x[ ,2]^2) / (50 * x[ ,1] * x[ ,2] + 10)
    return(out)
  }


###################################################
### code chunk number 47: laGP.Rnw:1653-1663
###################################################
library("tgp")
rect <- matrix(rep(0:1, 4), ncol = 2, byrow = 2)
ny <- 50
X <- lhs(ny, rect[1:2,] )
u <- c(0.2, 0.1)
Zu <- M(X, matrix(u, nrow = 1)) 
sd <- 0.5
reps <- 2
Y <- rep(Zu, reps) + rep(bias(X), reps) + 
  rnorm(reps * length(Zu), sd = sd) 


###################################################
### code chunk number 48: laGP.Rnw:1675-1685
###################################################
nz <- 10000
XU <- lhs(nz, rect)
XU2 <- matrix(NA, nrow=10 * ny, ncol = 4)
for(i in 1:10) {
  I <- ((i - 1) * ny + 1):(ny * i)
  XU2[I, 1:2] <- X
}
XU2[ ,3:4] <- lhs(10 * ny, rect[3:4, ])
XU <- rbind(XU, XU2)
Z <- M(XU[ ,1:2], XU[ ,3:4])


###################################################
### code chunk number 49: laGP.Rnw:1697-1701
###################################################
bias.est <- TRUE
methods <- rep("alc", 2)
da <- d <- darg(NULL, XU)
g <- garg(list(mle = TRUE), Y) 


###################################################
### code chunk number 50: laGP.Rnw:1712-1721
###################################################
beta.prior <- function(u, a = 2, b = 2, log = TRUE)
{
  if(length(a) == 1) a <- rep(a, length(u))
  else if(length(a) != length(u)) stop("length(a) must be 1 or length(u)")
  if(length(b) == 1) b <- rep(b, length(u))
  else if(length(b) != length(u)) stop("length(b) must be 1 or length(u)")
  if(log) return(sum(dbeta(u, a, b, log=TRUE)))
  else return(prod(dbeta(u, a, b, log=FALSE)))
}


###################################################
### code chunk number 51: laGP.Rnw:1729-1741
###################################################
initsize <- 10*ncol(X)
imesh <- 0.1
irect <- rect[1:2,]
irect[,1] <- irect[,1] + imesh/2
irect[,2] <- irect[,2] - imesh/2
uinit.cand <- lhs(10 * initsize, irect) 
uinit <- dopt.gp(initsize, Xcand = lhs(10 * initsize, irect))$XX
llinit <- rep(NA, nrow(uinit))
for(i in 1:nrow(uinit)) {
  llinit[i] <- fcalib(uinit[i,], XU, Z, X, Y, da, d, g, beta.prior, 
                  methods, M, bias.est, nth, verb = 0)
}


###################################################
### code chunk number 52: laGP.Rnw:1760-1763
###################################################
library("crs")
opts <- list("MAX_BB_EVAL" = 1000, "INITIAL_MESH_SIZE" = imesh, 
  "MIN_POLL_SIZE" = "r0.001", "DISPLAY_DEGREE" = 0)


###################################################
### code chunk number 53: laGP.Rnw:1777-1792
###################################################
its <- 0
o <- order(llinit)
i <- 1
out <- NULL
while(its < 10) {
  outi <- snomadr(fcalib, 2, c(0,0), 0, x0 = uinit[o[i],],
            lb = c(0,0), ub = c(1,1), opts = opts, XU = XU, 
            Z = Z, X = X, Y = Y, da = da, d = d, g = g, 
            methods = methods, M = M, bias = bias.est, 
            omp.threads = nth, uprior = beta.prior, 
            save.global = .GlobalEnv, verb = 0)
  its <- its + outi$iterations
  if(is.null(out) || outi$objective < out$objective) out <- outi
  i <- i + 1;
}


###################################################
### code chunk number 54: laGP.Rnw:1804-1809
###################################################
Xp <- rbind(uinit, as.matrix(fcalib.save[ ,1:2]))
Zp <- c(-llinit, fcalib.save[ ,3])
wi <- which(!is.finite(Zp))
if(length(wi) > 0) { Xp <- Xp[-wi, ]; Zp <- Zp[-wi]}
surf <- interp(Xp[ ,1], Xp[ ,2], Zp, duplicate = "mean")


###################################################
### code chunk number 55: usurf
###################################################
image(surf, xlab = "u1", ylab = "u2", main = "posterior surface",
  col = heat.colors(128), xlim = c(0,1), ylim = c(0,1))
points(uinit)
points(fcalib.save[,1:2], col = 3, pch = 18)
u.hat <- outi$solution
points(u.hat[1], u.hat[2], col = 4, pch = 18)
abline(v = u[2], lty = 2)
abline(h = u[1], lty = 2)


###################################################
### code chunk number 56: laGP.Rnw:1843-1848
###################################################
Xu <- cbind(X, matrix(rep(u, ny), ncol = 2, byrow = TRUE))
Mhat.u <- aGP.seq(XU, Z, Xu, da, methods, ncalib = 2, omp.threads = nth, 
  verb = 0)
cmle.u <- discrep.est(X, Y, Mhat.u$mean, d, g, bias.est, FALSE)
cmle.u$ll <- cmle.u$ll + beta.prior(u)


###################################################
### code chunk number 57: laGP.Rnw:1851-1852
###################################################
data.frame(u.hat = -outi$objective, u = cmle.u$ll)


###################################################
### code chunk number 58: laGP.Rnw:1864-1868
###################################################
nny <- 1000  
XX <- lhs(nny, rect[1:2,],)
ZZu <- M(XX, matrix(u, nrow = 1)) 
YYtrue <- ZZu + bias(XX) 


###################################################
### code chunk number 59: laGP.Rnw:1872-1879
###################################################
XXu <- cbind(XX, matrix(rep(u, nny), ncol = 2, byrow = TRUE))
Mhat.oos.u <- aGP.seq(XU, Z, XXu, da, methods, ncalib = 2, 
  omp.threads = nth, verb = 0)
YYm.pred.u <- predGP(cmle.u$gp, XX)
YY.pred.u <- YYm.pred.u$mean + Mhat.oos.u$mean
rmse.u <- sqrt(mean((YY.pred.u - YYtrue)^2))
deleteGP(cmle.u$gp)


###################################################
### code chunk number 60: laGP.Rnw:1886-1891
###################################################
Xu <- cbind(X, matrix(rep(u.hat, ny), ncol = 2, byrow = TRUE))
Mhat <- aGP.seq(XU, Z, Xu, da, methods, ncalib = 2, omp.threads = nth, 
  verb = 0)
cmle <- discrep.est(X, Y, Mhat$mean, d, g, bias.est, FALSE)
cmle$ll <- cmle$ll + beta.prior(u.hat)


###################################################
### code chunk number 61: laGP.Rnw:1896-1897
###################################################
print(c(cmle$ll, -outi$objective))


###################################################
### code chunk number 62: laGP.Rnw:1901-1907
###################################################
XXu <- cbind(XX, matrix(rep(u.hat, nny), ncol = 2, byrow = TRUE))
Mhat.oos <- aGP.seq(XU, Z, XXu, da, methods, ncalib = 2, 
  omp.threads = nth, verb = 0)
YYm.pred <- predGP(cmle$gp, XX)
YY.pred <- YYm.pred$mean + Mhat.oos$mean
rmse <- sqrt(mean((YY.pred - YYtrue)^2))


###################################################
### code chunk number 63: laGP.Rnw:1910-1911
###################################################
data.frame(u.hat = rmse, u = rmse.u)


