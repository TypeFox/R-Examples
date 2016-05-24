## Examples from Brockhoff, P.B. (2003). The statistical power of
## replications in difference tests.
## Food Quality and Preference, 14, pp. 405--417.
library(sensR)

## Data hunter 2, tab1, dataset 1:
x <- c(0, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4,
5, 5, 5, 6, 6, 6, 7, 7, 7, 8, 9, 9)
X <- matrix(c(x,rep(12,24)),ncol=2,byrow=F)
summary(bbc <- betabin(X, corr = TRUE, method = "triangle"))

target <- c(0.1175015, 0.2146566)
stopifnot(all(abs(coef(bbc)[1:2] - target) < 1e-5))
## Methods:
logLik(bbc)
vcov(bbc)
AIC(bbc)
bbc

sb <- summary(betabin(X, corr = TRUE, method = "twoAFC"))
b <- as.vector(coef(sb)[, 1])
b.expected <- c(0, 1, .5, 0, 0)
stopifnot(isTRUE(
    all.equal(b, b.expected, tol=1e-4)
    ))

## summary(betabin(X, corr = FALSE, pGuess = 1/3,
##                 method = "alpha-beta"))
## summary(betabin(X, corr = TRUE, pGuess = 1/3,
##                 method = "alpha-beta"))
sb <- summary(betabin(X, corr = TRUE, method = "threeAFC"))
print(coef(sb), digits=3)

## Hunter experiment 3, data set 2.
x <- c(2, 2, 3, 4, 5, 5, 5, 5, 6, 6, 6, 7,
       7, 7, 7, 7, 7, 8, 8, 8, 9, 11, 12)
X2 <- cbind(x, 12)
summary(betabin(X2, method = "triangle"))
summary(betabin(X2, cor = 0, method = "triangle"))
## summary(betabin(X2, pGuess=1/3, method = "a"))
## summary(betabin(X2, cor = 0, pGuess=1/3, method = "a"))
## Correspond to the estimates in Brockhoff(2003)

## Data set no. 3:
x <- c(33, 35, 36, 36, 98, 99)
X3 <- cbind(x, 100)
summary(betabin(X3, method = "triangle"))
summary(betabin(X3, cor = 0, method = "triangle"))
## summary(betabin(X3, pGuess=1/3, method = "a"))
## summary(betabin(X3, cor = 0, pGuess=1/3, method = "a"))
## The uncorrected estimates correspond to those in Brockhoff(2003),
## but the corrected ones do not.

## Hard example at the boundary:
## Data set no. 4:
x <- c(0, 1,1,1, 2,2,2, 3,3,3,3,3)
X4 <- cbind(x, 4)
sum(x)/(12*4) # 5
b <- coef(summary(betabin(X4, method = "triangle")))
stopifnot(isTRUE(
    all.equal(as.vector(b)[1:2], c(.25, 0), tol=1e-2)
    ))
b <- coef(summary(betabin(X4, cor = 0, method = "triangle")))
stopifnot(isTRUE(
    all.equal(as.vector(b)[1:2], c(.5, 0), tol=1e-4)
    ))
## summary(betabin(X4, pGuess=1/3, method = "a"))
## summary(betabin(X4, cor = 0, pGuess=1/3, method = "a"))

x <- c(0,0,1,1, 2,2, 3,3,3,3,3,3)
(X6 <- cbind(x, 4))
sum(x)/(12*4) # 5
summary(betabin(X6, method = "triangle"))
summary(betabin(X6, cor = 0, method = "triangle"))
## summary(betabin(X6, pGuess=1/3, method = "a"))
## summary(betabin(X6, cor = 0, pGuess=1/3, method = "a"))
## Correspond to the estimates in Brockhoff(2003)

## Data set no. 7:
x <- c(0,0,1,1, 2,2,2,2, 3,3,4,4)
(X7 <- cbind(x, 4))
sum(x)/(12*4) # 5
summary(betabin(X7, method = "triangle"))
summary(betabin(X7, cor = 0, method = "triangle"))
## summary(betabin(X7, pGuess=1/3, method = "a"))
## summary(betabin(X7, cor = 0, pGuess=1/3, method = "a"))
## The uncorrected estimates correspond to those in Brockhoff(2003),
## but the corrected ones do not.


## Data set no. 8:
x <- c(0,0,0,1, 2,2,2, 3,3,3,4,4)
(X8 <- cbind(x, 4))
sum(x)/(12*4) # 5
summary(betabin(X8, method = "triangle"))
summary(betabin(X8, method = "triangle", cor = 0))
## summary(betabin(X8, pGuess=1/3, method = "a"))
## summary(betabin(X8, cor = 0, pGuess=1/3, method = "a"))
## The uncorrected estimates correspond to those in Brockhoff(2003),
## but the corrected ones do not.


## Data set no. 9:
x <- c(0,0,0,1,1, 2,2, 3,3,4,4,4)
(X9 <- cbind(x, 4))
sum(x)/(12*4) # 5
summary(betabin(X9, method = "triangle"))
summary(betabin(X9, cor = 0, method = "triangle"))
## summary(betabin(X9, pGuess=1/3, method = "a"))
## summary(betabin(X9, cor = 0, pGuess=1/3, method = "a"))
## The uncorrected estimates correspond to those in Brockhoff(2003),
## but the corrected ones do not.

#################################
## Testing data sets from Bi(2006):

## Table 7.3, 2AFC, n=2, k=30
x <- c(0, 0, rep(1, 6), rep(2, 30-8))
X <- cbind(x, 2); X
summary(bbc <- betabin(X))
## Moment estimates: (mu, gamma) = (0.667, 0.700) (Bi, 2006)
## These ML estimates correspond to the ML
## estimates reported by Bi.
bbc$vcov
## The variance covariance matrix differs from the one reported by Bi,
## however.

## Table 7.4, 3AFC, n=4, k=30
x <- c(0, rep(2, 4), rep(3, 6), rep(4, 30-11))
X <- cbind(x, 4); X
summary(bbc <- betabin(X, method = "triangle"))
## Moment estimates: (mu, gamma) = (0.775, 0.4265) (Bi, 2006)
## ML and moment estimates are not far from each other
## These ML estimates correspond exactly to those by Bi(2006)
vcov(bbc)
## The variance-covariance matrix is not exactly the same however.
