## Particle Learning for regression Gaussian Processes on
## a simple 1-d sinusoidal function

## load the plgp library
library(plgp)
library(tgp) ## for the mean0.range1 function
## if mean0.range1 is not in the tgp NAMESPACE you may have
## to copy the source over manually

## close down old graphics windows and clear session
graphics.off()
rm(list=ls())

## set up 1-d Sinusoidal data following Gramacy & Lee (2007)
fsin <- function(X){ return((sin(pi*X/5) + 0.2*cos(4*pi*X/5))) }
rect <- c(0,9.6)
X <- lhs(50,rect)
Y <- fsin(X) + rnorm(length(X), sd=0.1)

## coerse the input matrix and adjust by the rectangle
## and put data into data-generating function
Xs <- rectscale(X, rect)
Ys <- mean0.range1(Y)$X  ## scale the response data
formals(data.GP)$X <- Xs
formals(data.GP)$Y <- Ys

## default prior
prior <- prior.GP(1, "isotropic")

## set up starting and ending times:
## by sepcifying start = end we can do Metropolis-Hastings
## for comparisons
start <- ncol(Xs) + 3
end <- nrow(Xs)

## Particle Learning Inference!
out <- PL(dstream=data.GP, ## static PL
          start=start, end=end,
          init=draw.GP,  ## initializing with Metropolis-Hastings
          lpredprob=lpredprob.GP, propagate=propagate.GP,
          prior=prior, addpall=addpall.GP,
          params=params.GP)

## design a grid of predictive locations
XX <- seq(rect[1], rect[2],length=99)
XXs <- rectscale(XX, rect)
## truth for RMSE comparisons
YY <- fsin(XX) 
YYs <- mean0.range1(YY)$X  

## sample from the particle posterior predictive distribution 
outp <- papply(XX=XXs, fun=pred.GP, Y=PL.env$pall$Y, quants=TRUE, prior=prior)

## unscale the data locations
X <- rectunscale(PL.env$pall$X, rect)

## set up to make two plots
par(mfrow=c(1,2))

## plot the individual lines of the predictive distribution
plot(X,Ys, main="predictive: each particle", xlab="x", ylab="y")
for(i in 1:length(outp)) {
  lines(XX, outp[[i]]$m)
  lines(XX, outp[[i]]$q1, col=2, lty=2)
  lines(XX, outp[[i]]$q2, col=2, lty=2)
}

## mean of the mean of the particles predictive
m <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) m <- m + outp[[i]]$m
m <- m / length(outp)

## a calculation of RMSE to the truth
rmse <- sqrt(mean(m - YYs)^2)
rmse

## mean of the quantiles of the particles predictive
q2 <- q1 <- rep(0, nrow(as.matrix(XX)))
for(i in 1:length(outp)) {
  q1 <- q1 + outp[[i]]$q1
  q2 <- q2 + outp[[i]]$q2
}
q1 <- q1 / length(outp)
q2 <- q2 / length(outp)

## plot the summary stats of the predictive distribution
plot(X, Ys, main="predictive surface",
     ylim=range(c(q1,q2)), xlab="x", ylab="y")
lines(XX, m, lwd=2)
lines(XX, q1, col=2, lty=2, lwd=2)
lines(XX, q2, col=2, lty=2, lwd=2)
legend("topright", c("mean", "90% quantiles"), lty=1:2, col=1:2, lwd=2)

## look at a historgram of the parameters
params <- params.GP()
dev.new()
par(mfrow=c(1,2)) ## two plots
hist(params$d)
hist(params$g)

## make a plot of the joint samples of d and g
par(mfrow=c(1,1))  ## back to one plot
plot(params$d, params$g, xlab="d", ylab="g",
     main="Samples of range (d) and nugget (g)")
