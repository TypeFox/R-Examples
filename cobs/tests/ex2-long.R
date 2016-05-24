####
suppressMessages(library(cobs))

options(digits = 5)
if(!dev.interactive(orNone=TRUE)) pdf("ex2.pdf")

source(system.file("util.R", package = "cobs"))

set.seed(821)
x <- round(sort(rnorm(200)), 3) # rounding -> multiple values
sum(duplicated(x)) # 9
y <- (fx <- exp(-x)) + rt(200,4)/4
summaryCobs(cxy  <- cobs(x,y, "decrease"))
1 - sum(cxy $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 97.6%

## Interpolation
cpuTime(cxyI  <- cobs(x,y, "decrease", knots = unique(x)))
## takes quite long : 63 sec. (Pent. III, 700 MHz) --- this is because
## each knot is added sequentially...  {{improve!}}
summaryCobs(cxyI)# only 7 knots remaining!

summaryCobs(cxy1  <- cobs(x,y, "decrease", lambda = 0.1))
1 - sum(cxy1 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 98.2%

summaryCobs(cxy2 <- cobs(x,y, "decrease", lambda = 1e-2))
1 - sum(cxy2 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 98.2% (tiny bit better)

summaryCobs(cxy3 <- cobs(x,y, "decrease", lambda = 1e-6, nknots = 60))
1 - sum(cxy3 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 98.36%

cpuTime(cxy4 <- cobs(x,y, "decrease", lambda = 1e-6, nknots = 100))# ~ 3 sec.
1 - sum(cxy4 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 98.443%

cpuTime(cxy5 <- cobs(x,y, "decrease", lambda = 1e-6, nknots = 150))# ~ 8.7 sec.
1 - sum(cxy5 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 98.4396%

## regularly spaced x :
X <- seq(-1,1, len = 201)
xx <- c(seq(-1.1, -1, len = 11), X,
        seq( 1,  1.1, len = 11))
y <- (fx <- exp(-X)) + rt(201,4)/4
summaryCobs(cXy  <- cobs(X,y, "decrease"))
1 - sum(cXy $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 77.2%

(cXy.9 <- cobs(X,y, "decrease", tau = 0.9))
(cXy.1 <- cobs(X,y, "decrease", tau = 0.1))
(cXy.99<- cobs(X,y, "decrease", tau = 0.99))
(cXy.01<- cobs(X,y, "decrease", tau = 0.01))
plot(X,y, xlim = range(xx),
     main = "cobs(*, \"decrease\"), N=201, tau = 50% (Med.), 1,10, 90,99%")
lines(predict(cXy,    xx), col = 2)
lines(predict(cXy.1,  xx), col = 3)
lines(predict(cXy.9,  xx), col = 3)
lines(predict(cXy.01, xx), col = 4)
lines(predict(cXy.99, xx), col = 4)

## Interpolation
cpuTime(cXyI  <- cobs(X,y, "decrease", knots = unique(X)))
## takes ~ 47 sec. (Pent. III, 700 MHz)
summaryCobs(cXyI)# only 7 knots remaining!

summaryCobs(cXy1 <- cobs(X,y, "decrease", lambda= 0.1))
1 - sum(cXy1 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 79.53 %

summaryCobs(cXy2 <- cobs(X,y, "decrease", lambda= 1e-2))
1 - sum(cXy2 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 80.004%

summaryCobs(cXy3 <- cobs(X,y, "decrease", lambda= 1e-6, nknots = 60))
1 - sum(cXy3 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 80.424%

cpuTime(cXy4 <- cobs(X,y, "decrease", lambda= 1e-6, nknots = 100))#~16.5"
## not converged (in 4020 iter.)
1 - sum(cXy4 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 80.517%

cpuTime(cXy5 <- cobs(X,y, "decrease", lambda= 1e-6, nknots = 150))#~12.8"
1 - sum(cXy5 $ resid ^ 2) / sum((y - mean(y))^2) # R^2 = 81.329%







