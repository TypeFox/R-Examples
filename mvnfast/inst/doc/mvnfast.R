## ----setup, include=FALSE------------------------------------------------
library(knitr)
opts_chunk$set(out.extra='style="display:block; margin: auto"', fig.align="center", tidy=FALSE)

## ----rmvn----------------------------------------------------------------
library("microbenchmark")
library("mvtnorm")
library("mvnfast")
library("MASS")

N <- 10000
d <- 20

# Creating mean and covariance matrix
mu <- 1:d
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)

microbenchmark(rmvn(N, mu, mcov, ncores = 2),
               rmvn(N, mu, mcov),
               rmvnorm(N, mu, mcov),
               mvrnorm(N, mu, mcov))

## ----dmvn----------------------------------------------------------------
# Generating random vectors 
N <- 10000
d <- 20
mu <- 1:d
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)
X <- rmvn(N, mu, mcov)

microbenchmark(dmvn(X, mu, mcov, ncores = 2),
               dmvn(X, mu, mcov),
               dmvnorm(X, mu, mcov))

## ----maha----------------------------------------------------------------
# Generating random vectors 
N <- 10000
d <- 20
mu <- 1:d
tmp <- matrix(rnorm(d^2), d, d)
mcov <- tcrossprod(tmp, tmp)
X <- rmvn(N, mu, mcov)

microbenchmark(maha(X, mu, mcov, ncores = 2),
               maha(X, mu, mcov),
               mahalanobis(X, mu, mcov))

## ----mixSim--------------------------------------------------------------
set.seed(5135)
N <- 10000
d <- 2
mu1 <- c(0, 0); mu2 <- c(2, 3)
Cov1 <- matrix(c(1, 0, 0, 2), 2, 2)
Cov2 <- matrix(c(1, -0.9, -0.9, 1), 2, 2)

bin <- rbinom(N, 1, 0.5)

X <- bin * rmvn(N, mu1, Cov1) + (!bin) * rmvn(N, mu2, Cov2)

## ----mixPlot-------------------------------------------------------------
# Plotting
np <- 100
xvals <- seq(min(X[ , 1]), max(X[ , 1]), length.out = np)
yvals <- seq(min(X[ , 2]), max(X[ , 2]), length.out = np)
theGrid <- expand.grid(xvals, yvals) 
theGrid <- as.matrix(theGrid)
dens <- 0.5 * dmvn(theGrid, mu1, Cov1) + 0.5 * dmvn(theGrid, mu2, Cov2)
plot(X[ , 1], X[ , 2], pch = '.', lwd = 0.01, col = 3)
contour(x = xvals, y = yvals, z = matrix(dens, np, np),
        levels = c(0.002, 0.01, 0.02, 0.04, 0.08, 0.15 ), add = TRUE, lwd = 2)

# Mean-shift
library(plyr)
inits <- matrix(c(-2, 2, 0, 3, 4, 3, 2, 5, 2, -3, 2, 2, 0, 2, 3, 0, 0, -4, -2, 6), 
                10, 2, byrow = TRUE)
traj <- alply(inits,
              1,
              function(input)
                  ms(X = X, 
                     init = input, 
                     H = 0.05 * cov(X), 
                     ncores = 2, 
                     store = TRUE)$traj
              )

invisible( lapply(traj, 
                  function(input){ 
                    lines(input[ , 1], input[ , 2], col = 2, lwd = 1.5)
                    points(tail(input[ , 1]), tail(input[ , 2]))
           }))

