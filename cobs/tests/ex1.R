#### OOps! Running this in 'CMD check' or in *R* __for the first time__
#### ===== gives a wrong result (at the end) than when run a 2nd time
####-- problem disappears with introduction of   if (psw) call ... in Fortran

suppressMessages(library(cobs))
options(digits = 6)
if(!dev.interactive(orNone=TRUE)) pdf("ex1.pdf")

source(system.file("util.R", package = "cobs"))

## Simple example from  example(cobs)
set.seed(908)
x <- seq(-1,1, len = 50)
f.true <- pnorm(2*x)
y <- f.true + rnorm(50)/10
## specify constraints (boundary conditions)
con <- rbind(c( 1,min(x),0),
             c(-1,max(x),1),
             c( 0, 0,  0.5))
## obtain the median *regression* B-spline using automatically selected knots
coR <- cobs(x,y,constraint = "increase", pointwise = con)
summaryCobs(coR)
coR1 <- cobs(x,y,constraint = "increase", pointwise = con, degree = 1)
summary(coR1)

## compute the median *smoothing* B-spline using automatically chosen lambda
coS <- cobs(x,y,constraint = "increase", pointwise = con,
            lambda = -1, trace = 3)
with(coS, cbind(pp.lambda, pp.sic, k0, ifl, icyc))
with(coS, plot(pp.sic ~ pp.lambda, type = "b", log = "x", col=2,
	       main = deparse(call)))
##-> very nice minimum close to 1

summaryCobs(coS)

plot(x, y, main = "cobs(x,y, constraint=\"increase\", pointwise = *)")
matlines(x, cbind(fitted(coR), fitted(coR1), fitted(coS)),
         col = 2:4, lty=1)

##-- real data example (still n = 50)
data(cars)
attach(cars)
co1   <- cobs(speed, dist, "increase")
co1.1 <- cobs(speed, dist, "increase", knots.add = TRUE)
co1.2 <- cobs(speed, dist, "increase", knots.add = TRUE, repeat.delete.add = TRUE)
## These three all give the same -- only remaining knots (outermost data):
ic <- which("call" == names(co1))
stopifnot(all.equal(co1[-ic], co1.1[-ic]),
          all.equal(co1[-ic], co1.2[-ic]))
1 - sum(co1   $ resid ^2) / sum((dist - mean(dist))^2) # R^2 = 64.2%

co2 <- cobs(speed, dist, "increase", lambda = -1)# 6 warnings
summaryCobs(co2)
1 - sum(co2 $ resid ^2) / sum((dist - mean(dist))^2)# R^2= 67.4%

co3 <- cobs(speed, dist, "convex", lambda = -1)# 3 warnings
summaryCobs(co3)
1 - sum(co3 $ resid ^2) / sum((dist - mean(dist))^2) # R^2 = 66.25%

with(co2, plot(pp.sic ~ pp.lambda, type = "b", col = 3, log = "x"))
with(co3, plot(pp.sic ~ pp.lambda, type = "b", col = 4, log = "x"))

plot(speed, dist, main = "cobs(speed,dist, ..) for data(cars)")
lines(speed, fitted(co2),   col=3); rug(knots(co2),   col=3)
lines(speed, fitted(co3),   col=4); rug(knots(co3),   col=4)
lines(speed, fitted(co1.1), col=2); rug(knots(co1.1), col=2)
detach(cars)

##-- another larger example using "random" x
set.seed(101)
x <- round(sort(rnorm(500)), 3) # rounding -> multiple values
sum(duplicated(x)) # 32
y <- (fx <- exp(-x)) + rt(500,4)/4
summaryCobs(cxy  <- cobs(x,y, "decrease"))
1 - sum(cxy $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 95.9%

## Interpolation
if(FALSE) { ##-- since it takes too long here!
   cpuTime(cxyI  <- cobs(x,y, "decrease", knots = unique(x)))
   ## takes very long : 1864.46 sec. (Pent. III, 700 MHz)
   summaryCobs(cxyI)# only 8 knots remaining!
}

dx <- diff(range(ux <- unique(x)))
rx <- range(xx <- seq(ux[1] - dx/20, ux[length(ux)] + dx/20, len = 201))
cpuTime(cxyI  <- cobs(x,y, "decrease", knots = ux, nknots = length(ux)))
## 17.3 sec. (Pent. III, 700 MHz)
summary(cxyI)
pxx <- predict(cxyI, xx)
plot(x,y, cex = 3/4, xlim = rx, ylim = range(y, pxx[,"fit"]),
     main = "Artificial (x,y), N=500 : `interpolating' cobs()")
lines(xx, exp(-xx), type = "l", col = "gray40")
lines(pxx, col = "red")
rug(cxyI$knots, col = "blue", lwd = 0.5)

## Deg = 1
cpuTime(cI1 <- cobs(x,y, "decrease", knots= ux, nknots= length(ux), degree = 1))
summary(cI1)
pxx <- predict(cI1, xx)
plot(x,y, cex = 3/4, xlim = rx, ylim = range(y, pxx[,"fit"]),
     main = paste("Artificial, N=500, `interpolate'", deparse(cI1$call)))
lines(xx, exp(-xx), type = "l", col = "gray40")
lines(pxx, col = "red")
rug(cI1$knots, col = "blue", lwd = 0.5)


cpuTime(cxyS <- cobs(x,y, "decrease", lambda = -1))
## somewhat <  2 sec. (Pent. III, 700 MHz)
pxx <- predict(cxyS, xx)
pxx[xx > max(x) , ]# those outside to the right -- currently all = Inf !
summaryCobs(cxyS)
R2 <- 1 - sum(cxyS $ resid ^ 2) / sum((y - mean(y))^2)
R2 # R^2 = 96.3%, now 96.83%

plot(x,y, cex = 3/4, xlim = rx, ylim = range(y, pxx[,"fit"], finite = TRUE),
     main = "Artificial (x,y), N=500 : cobs(*, lambda = -1)")
mtext(substitute(R^2 == r2 * "%", list(r2 = round(100*R2,1))))
lines(xx, exp(-xx), type = "l", col = "gray40")
lines(pxx, col = "red")
rug(cxyS$knots, col = "blue", lwd = 1.5)

## Show print-monitoring :

cxyS <- cobs(x,y, "decrease", lambda = -1, print.mesg = 2)# << improve! (1 line)
cxyS <- cobs(x,y, "none",     lambda = -1, print.mesg = 3)

## this does NOT converge (and "trace = 3" does *not* show it -- improve!)

cxyC <- cobs(x,y, "concave", lambda = -1)
summaryCobs(cxyC)


dev.off()
