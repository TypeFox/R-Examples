
library("mlt")
library("truncreg")


## Left-truncated 
set.seed(29)
n <- 1000	
x <- runif(n, max = 2 * pi)
y <- rnorm(n, 2*x + 1, .25)
d <- data.frame(y = y, x = x)
## truncated response
tr <- 2.5
d$yt <- ifelse(d$y > tr, d$y, NA)
d$trunc_left <- tr

tmod <- truncreg(yt ~ x, data = d, point = tr, direction = "left")
coef(tmod)
logLik(tmod)
vcov(tmod)

## MLT
fm <- as.formula("~ x")
yb <- polynomial_basis(numeric_var("yt", support = range(d$yt, na.rm = TRUE)), 
                       coef = c(TRUE, TRUE), ci = c(-Inf, 0))
m <- ctm(yb, shifting = fm, todistr = "Normal", data = d)
dfull <- d[!is.na(d$yt),,drop = FALSE]
dfull$yt <- R(dfull$yt, tleft = dfull$trunc_left)
mltmod <- mlt(m, data = dfull)
(cf <- coef(mltmod))
c(-cf[1] / cf[2], -cf[3] / cf[2], 1 / cf[2])
logLik(mltmod)
vcov(mltmod)

#library("numDeriv")

solve(numDeriv::hessian(mltmod$loglik, coef(mltmod), 
                        weights = weights(mltmod)))


## right-truncated

set.seed(29)
n <- 1000	
x <- runif(n, max = 2 * pi)
y <- rnorm(n, 2*x + 1, .25)
d <- data.frame(y = y, x = x)
## truncated response
tr <- 11
d$yt <- ifelse(d$y < tr, d$y, NA)
d$trunc_right <- tr

tmod <- truncreg(yt ~ x, data = d, point = tr, direction = "right")
coef(tmod)
logLik(tmod)
vcov(tmod)

## MLT
fm <- as.formula("~ x")
yb <- polynomial_basis(numeric_var("yt", support = range(d$yt, na.rm = TRUE)), 
                       coef = c(TRUE, TRUE), ci = c(-Inf, 0))
m <- ctm(yb, shifting = fm, todistr = "Normal", data = d)
dfull <- d[!is.na(d$yt),,drop = FALSE]
dfull$yt <- R(dfull$yt, tright = dfull$trunc_right)
mltmod <- mlt(m, data = dfull)
(cf <- coef(mltmod))
c(-cf[1] / cf[2], -cf[3] / cf[2], 1 / cf[2])
logLik(mltmod)
vcov(mltmod)

#library("numDeriv")

solve(numDeriv::hessian(mltmod$loglik, coef(mltmod), 
                        weights = weights(mltmod)))


