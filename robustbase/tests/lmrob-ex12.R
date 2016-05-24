
library(robustbase)
source(system.file("test-tools-1.R", package="Matrix", mustWork=TRUE))
##-> assertError(), etc

set.seed(1) # since now .Random.seed is used by default!

## EX 1
data(coleman)
## "empty model" (not really a lot of sense)
(m0 <- lmrob(Y ~ 0, data = coleman))
summary(m0)
## "Intercept" only: robust mean
(m1 <- lmrob(Y ~ 1, data = coleman))
summary(m1)


(mC <- lmrob(Y ~ ., data = coleman,
	     control = lmrob.control(refine.tol = 1e-8, rel.tol = 1e-9)))
summary(mC)
## Values will change once we use R's random number generator !
stopifnot(
all.equal(unname(coef(mC)),
	  c(30.50232, -1.666147, 0.08425381, 0.6677366, 1.167777, -4.136569),
	  tolerance = 2e-7)# 6.112 e-8 (32-b)
)
dput(signif(unname(coef(mC)), 7))
## 64b(0.2-0): c(30.50232, -1.666147, 0.08425381, 0.6677366, 1.167777, -4.136569)
## 32b(0.2-0):	 "exactly" same !
## Full precision:
dput(unname(coef(mC)))
## 2012-06-04:
## 32-bit:c(30.5023184450149, -1.66614687548007, 0.0842538074792178, 0.667736590070332, 1.16777744029117, -4.13656885405815)
## 64-bit:c(30.5023184450148, -1.66614687548008, 0.0842538074792178, 0.667736590070332, 1.16777744029117, -4.13656885405814)
##
## 32-bit:c(30.5023183940104, -1.66614687550933, 0.0842538074635567, 0.667736589938547, 1.16777744089398, -4.13656884777543)
## 64-bit:c(30.5023184150851, -1.66614687537736, 0.0842538074722959, 0.667736589980183, 1.16777744061092, -4.1365688503035)

str(mC)

## EX 2
gen <- function(n,p, n0, y0, x0, beta = rep(1, p))
{
    stopifnot(n >= 1, p >= 1, n0 >= 0, length(beta) == p)
    x <- matrix(rnorm(n*p),n,p) # iid x's
    y <- x %*% beta + rnorm(n)
    xc <- matrix(0,n0,p)
    xc[,1] <- x0
    xc <- xc + 0.1*matrix(rnorm(n0*p),n0,p)
    x[1:n0,] <- xc
    y[1:n0] <- y0 + .001*rnorm(n0)
    list(x=x, y=y)
}

## generate --a sample of  n  observations with	 p  variables
## and 10% of outliers near (x1,y) = (10,10)
n <- 500 ; n0 <- n %/% 10
p <- 7 ## p = 20 is more impressive but too slow for "standard test"
set.seed(17)
a <- gen(n=n, p=p, n0= n0, y0=10, x0=10)
plot(a$x[,1], a$y, col = c(rep(2, n0), rep(1, n-n0)))
system.time( m1 <- lmrob(y~x, data = a,
                         control = lmrob.control(compute.rd = TRUE)))
plot(m1, ask=FALSE)
##-> currently 5 plots; MM:I  don't like #3 (Response vs fitted)

## don't compute robust distances --> faster by factor of two:
system.time(m2 <- lmrob(y~x, data = a,
			control = lmrob.control(compute.rd = FALSE)))
## ==> half of the CPU time is spent in covMcd()!
(sm2 <- summary(m2))
l1 <- lm(y~x, data = a)
cbind(robust = coef(sm2)[,1:2],
      lm = coef(summary(l1))[,1:2])

m2.S1 <- with(a, lmrob.S(cbind(1,x), y, trace.lev = 2,
			 ## trace.lev = 2 : quite a bit of output
			 control= lmrob.control(seed = .Random.seed,
			 nRes = 80, k.max = 20, refine.tol = 1e-4)))
S.ctrl <- lmrob.control(seed = .Random.seed,## << keeps .Random.seed unchanged
			nResample = 1000, best.r.s = 15, refine.tol = 1e-9)
m2.S <- with(a, lmrob.S(cbind(1,x), y, control = S.ctrl, trace.lev = 1))
str(m2.S)

##--- Now use n > 2000 --> so we use C internal fast_s_large_n(...)
n <- 2500 ; n0 <- n %/% 10
a2 <- gen(n=n, p = 3, n0= n0, y0=10, x0=10)
plot(a2$x[,1], a2$y, col = c(rep(2, n0), rep(1, n-n0)))
rs <- .Random.seed
system.time( m3 <- lmrob(y~x, data = a2) )
m3
nrs <- .Random.seed # <-- to check that using 'seed' keeps .Random.seed
system.time( m4 <- lmrob(y~x, data = a2, seed = rs, compute.rd = FALSE))
(sm4 <- summary(m4))

## random seed must be the same because we used	 'seed = *' :
stopifnot(nrs == .Random.seed, identical(coef(m3), coef(m4)))

dput(signif(cf <- unname(coef(m3)), 7))
## 2012-06-04:c(-0.05108914, 1.005971, 1.003201, 0.9833263) - 32 AND 64 bit
##
## 0.2-0: c(0.007446546, 1.000712, 1.027921, 0.9896527)
## 0.2-1: c(0.03148659, 0.9980933, 1.016364, 1.03243)
## both for 32 and 64 bit

dput(signif(100 * (sd <- unname(coef(sm4)[, "Std. Error"])), 7))
## 2012-06-04:c(2.213815, 0.2864678, 2.202318, 2.180886) - 32 AND 64 bit
##
## 0.2-0: c(2.219388, 0.274644,  2.196982, 2.26253)
## 0.2-1: c(2.194914, 0.2737579, 2.371728, 2.206261)
## both for 32 and 64 bit

stopifnot(
	  all.equal(cf, c(-0.05108914, 1.00597115, 1.00320052, 0.98332632), tolerance= 7e-7)
	  , # ... e-7	 needed on 64b
	  all.equal(100*sd,c(2.2138147, 0.2864678, 2.2023182, 2.1808862),tolerance= 7e-7)
	  ) # 1.334 e-7	 needed on 64b

cat('Time elapsed: ', proc.time(),'\n') # "stats"

## rm(a,m1, m2, m3, m4, sm2, l1)

## Small examples from R-SIG-robust

## First example from René Locher :
dat1 <- data.frame(lconc= log(c(21.8, 23.7, 12.2, 38.5, 21, 38.9)),
                   dist =     c( 100, 180,  280,  30,  220,  6))
m5 <- lmrob(lconc ~ dist, data = dat1)
## Warning messages:
## ... S refinements did not converge (to tol=1e-07) in 200 iterations
##  "  "     "
m5$init.S$converged # FALSE
m5. <- lmrob(lconc ~ dist, data = dat1,
             control = lmrob.control(refine.tol = 1e-5))
m5.$init.S$converged # TRUE
## gives TRUE as the IRWLS iterations after the lmrob.S() have converged.

## 2nd example from René Locher , 6 Jun 2007

dat2 <- data.frame(lconc=log(c(29.5,40.1,21.1,25.3,27.3,25.2,26.9,19.1,16.4)),
                     dist =    c(520, 1480,1780, 740, 540,1050,1100,1640,1860))

res2 <- lmrob(lconc~dist, data = dat2)

## Used to give Warning messages:
## 1: rwls(): not converged in 1000 lambda iterations
## ...
## 4: rwls(): ............

res2 <- lmrob(lconc~dist, data = dat2, trace.lev = 3)
##                                     -------------
summary(res2)
stopifnot(dim(model.matrix(res2)) == c(9,2))

## Check predict():
dd <- seq(300, 2000, by = 50)
with(dat2, plot(dist, lconc, pch=20, cex=2, xlim = range(dd)))
new.d <- data.frame(dist=dd)
fit.dd <- predict(res2, new.d)
lines(dd, fit.dd, col=2, type="o")
predict(res2, new.d, se=TRUE)$se.fit
matlines(dd, predict(res2, new.d, interval="confidence")[, 2:3], col=3)

## Check handling of X of not full rank
test <- function(n, ...) {
    X <- matrix(c(rep(1:3, length.out = n), rnorm(2*n)), n, 4)
    y <- rnorm(n)
    X[,4] <- X[,2] + X[,3]
    X <- data.frame(X)
    X$X1 <- factor(X$X1)
    fail <- suppressWarnings(try(lmrob(y ~ ., X, ...), silent=TRUE))
    stopifnot(is(fail, "lmrob"))
}
set.seed(0)
test(12) ## fast_S()
test(2500) ## fast_S_large_n()
test(200, trace.lev = TRUE)

## Check a case, where cov() matrix needs "posdefify":

coleman16 <- coleman[ -c(2, 7, 16, 19),]
(m16 <- lmrob(Y ~ ., data = coleman16, tuning.psi = 3.44, trace.lev = TRUE))
## failed in 0.9_0

assertWarning(
 lmrob(Y ~ ., data = coleman, setting = "KS2011", control = lmrob.control())
)

cat('Time elapsed: ', proc.time(),'\n') # "stats"
