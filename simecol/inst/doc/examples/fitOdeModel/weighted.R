library(simecol)

## ======== load example model =========
data(chemostat)

#source("chemostat.R")

## derive scenario
cs1 <- chemostat

## generate some noisy data
parms(cs1)[c("vm", "km")] <- c(2, 10)
times(cs1) <- c(from=0, to=20, by=2)
yobs <- out(sim(cs1))
obstime <- yobs$time
yobs$time <- NULL
yobs$S <- yobs$S + rnorm(yobs$S, sd= 0.1 * sd(yobs$S))*2
yobs$X <- yobs$X + rnorm(yobs$X, sd= 0.1 * sd(yobs$X))

## ======== optimize it! =========

times(cs1) <- obstime
solver(cs1) <- "lsoda"

whichpar  <- c("vm", "km")
parms(cs1)[whichpar] <- c(vm=1, km=2)

lower <- c(vm=0, km=0)
upper <- c(vm=4, km=20)

yobs$X[6]  <- 100 # outlier
weights    <- data.frame(X=rep(1, nrow(yobs)), S=rep(1, nrow(yobs)))
weights$X[6] <- .1

## non-weighted
res <- fitOdeModel(cs1, whichpar=c("vm", "km"), obstime, yobs,
  debuglevel=0, fn = ssqOdeModel,
  method = "PORT", lower = lower, upper = upper,
  #weights = weights,
  control=list(trace=TRUE),
  atol=1e-4, rtol=1e-4)

## assign fitted parameters to scenario cs1
parms(cs1)[whichpar] <- res$par

## weighted
cs2 <- cs1 # copy of the simulation object

res <- fitOdeModel(cs2, whichpar=c("vm", "km"), obstime, yobs,
  debuglevel=0, fn = ssqOdeModel,
  method = "PORT", lower = lower, upper = upper,
  weights = weights,
  control=list(trace=TRUE),
  atol=1e-4, rtol=1e-4)

## assign fitted parameters to scenario cs2
parms(cs2)[whichpar] <- res$par

## set small external time step for good graphics and simulate again
times(cs1) <- c(from=0, to=20, by=.1)
times(cs2) <- c(from=0, to=20, by=.1)
ysim1 <- out(sim(cs1))
ysim2 <- out(sim(cs2))

## compare results
par(mfrow=c(2,1))
plot(obstime, yobs$X, ylim = range(yobs$X, ysim1$X))
lines(ysim1$time, ysim1$X, col="blue")
lines(ysim2$time, ysim2$X, col="red")
plot(obstime, yobs$S, ylim= range(yobs$S, ysim1$S))
lines(ysim1$time, ysim1$S, col="blue")
lines(ysim2$time, ysim2$S, col="red")
