suppressMessages(library(cobs))

source(system.file("util.R", package = "cobs"))

data(DublinWind)
attach(DublinWind)##-> speed & day (instead of "wind.x" & "DUB.")
iday <- sort.list(day)

if(!dev.interactive(orNone=TRUE)) pdf("wind.pdf", width=10)

stopifnot(identical(day,c(rep(c(rep(1:365,3),1:366),4),
                          rep(1:365,2))))
co50.1 <- cobs(day, speed, constraint= "periodic", tau= .5, lambda= 2.2,
               degree = 1)
co50.2 <- cobs(day, speed, constraint= "periodic", tau= .5, lambda= 2.2,
               degree = 2)

plot(day,speed)
lines(day[iday], fitted(co50.1)[iday], col="orange", lwd = 2)
lines(day[iday], fitted(co50.2)[iday], col="light blue", lwd = 2)
rug(knots(co50.1), col=3)

nknots <- 13
## Compute the quadratic median smoothing B-spline using SIC
## lambda selection
co.o50 <-
 cobs(day, speed, knots.add = TRUE, constraint="periodic", nknots = nknots,
      tau = .5, lambda = -1, factor = 3, method = "uniform")
summary(co.o50)

op <- par(mfrow = c(3,1), mgp = c(1.5, 0.6,0), mar=.1 + c(3,3:1))
with(co.o50, plot(pp.sic ~ pp.lambda, type ="o",
                  col=2, log = "x", main = "co.o50: periodic"))
with(co.o50, plot(pp.sic ~ pp.lambda, type ="o", ylim = robrng(pp.sic),
                  col=2, log = "x", main = "co.o50: periodic"))
of <- 0.64430538125795
with(co.o50, plot(pp.sic - of ~ pp.lambda, type ="o", ylim = c(6e-15, 8e-15),
                  ylab = paste("sic -",formatC(of, dig=14, small.m = "'")),
                  col=2, log = "x", main = "co.o50: periodic"))
par(op)

## cobs99: Since SIC chooses a lambda that corresponds to the smoothest
## possible fit, rerun cobs with a larger lstart value
## (lstart <- log(.Machine$double.xmax)^3) # 3.57 e9
##
co.o50. <-
    cobs(day,speed,knots.add = TRUE, constraint = "periodic", nknots = 10,
         tau = .5,lambda = -1, method = "quantile")
summary(co.o50.)
(pc.5 <- predict(co.o50., interval = "both"))

co.o50.. <- cobs(day,speed, knots.add = TRUE, repeat.delete.add=TRUE,
                 constraint = "periodic", nknots = 10,
                 tau = .5, lambda = -1, method = "quantile")

co.o9 <- ## Compute the .9 quantile smoothing B-spline
    cobs(day,speed,knots.add = TRUE, constraint = "periodic", nknots = 10,
         tau = .9,lambda = -1, method = "uniform")
summary(co.o9)
(pc.9 <- predict(co.o9,interval = "both"))

co.o1 <- ## Compute the .1 quantile smoothing B-spline
    cobs(day,speed,knots.add = TRUE, constraint = "periodic",nknots = nknots,
         tau = .1,lambda = -1, method = "uniform")
summary(co.o1)
(pc.1 <- predict(co.o1, interval = "both"))

par(mfrow = c(1,2), mgp = c(1.5, .6,0), mar = .1 + c(3,3,1,1))
plot(day,speed, pch = 3, cex=0.6,  xlab = "DAYS", ylab = "SPEED (knots)")
lines(pc.5, lwd = 2.5, col = 2)
lines(pc.9, lwd = 2., col = "blue")
lines(pc.1, lwd = 2., col = "blue")
plot(day,speed,type = "n",xlab = "DAYS", ylab = "SPEED (knots)")
lines(pc.5, lwd = 1.5)
lines(pc.9, col = 3)
lines(pc.1, col = 3)
abline(v = co.o50.$knots, lty = 3, col = "gray70")
## rather rug(co.o5$knots, lty = 2)
