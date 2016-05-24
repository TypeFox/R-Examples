# pkgname <- "FKF"
# source(file.path(R.home("share"), "R", "examples-header.R"))
# options(warn = 1)
# library('FKF')
# 
# assign(".oldSearch", search(), pos = 'CheckExEnv')
# cleanEx()
# nameEx("fkf")
# ### * fkf
# 
# flush(stderr()); flush(stdout())

### Name: fkf
### Title: Fast Kalman filter
### Aliases: fkf
### Keywords: algebra models multivariate

### ** Examples

## <--------------------------------------------------------------------------->
## Example 1: ARMA(2, 1) model estimation.
## <--------------------------------------------------------------------------->
## This example shows how to fit an ARMA(2, 1) model using this Kalman
## filter implementation (see also stats' makeARIMA and KalmanRun).
n <- 1000

## Set the AR parameters
ar1 <- 0.6
ar2 <- 0.2
ma1 <- -0.2
sigma <- sqrt(0.2)

## Sample from an ARMA(2, 1) process
a <- arima.sim(model = list(ar = c(ar1, ar2), ma = ma1), n = n,
               innov = rnorm(n) * sigma)

## Create a state space representation out of the four ARMA parameters
arma21ss <- function(ar1, ar2, ma1, sigma) {
    Tt <- matrix(c(ar1, ar2, 1, 0), ncol = 2)
    Zt <- matrix(c(1, 0), ncol = 2)
    ct <- matrix(0)
    dt <- matrix(0, nrow = 2)
    GGt <- matrix(0)
    H <- matrix(c(1, ma1), nrow = 2) * sigma
    HHt <- H %*% t(H)
    a0 <- c(0, 0)
    P0 <- matrix(1e6, nrow = 2, ncol = 2)
    return(list(a0 = a0, P0 = P0, ct = ct, dt = dt, Zt = Zt, Tt = Tt, GGt = GGt,
                HHt = HHt))
}

## The objective function passed to 'optim'
objective <- function(theta, yt) {
    sp <- arma21ss(theta["ar1"], theta["ar2"], theta["ma1"], theta["sigma"])
    ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
               Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = yt)
    return(-ans$logLik)
}

theta <- c(ar = c(0, 0), ma1 = 0, sigma = 1)
fit <- optim(theta, objective, yt = rbind(a), hessian = TRUE)
fit

## Confidence intervals
rbind(fit$par - qnorm(0.975) * sqrt(diag(solve(fit$hessian))),
      fit$par + qnorm(0.975) * sqrt(diag(solve(fit$hessian))))

## Filter the series with estimated parameter values
sp <- arma21ss(fit$par["ar1"], fit$par["ar2"], fit$par["ma1"], fit$par["sigma"])
ans <- fkf(a0 = sp$a0, P0 = sp$P0, dt = sp$dt, ct = sp$ct, Tt = sp$Tt,
           Zt = sp$Zt, HHt = sp$HHt, GGt = sp$GGt, yt = rbind(a))

## Compare the prediction with the realization
plot(ans, at.idx = 1, att.idx = NA, CI = NA)
lines(a, lty = "dotted")

## Compare the filtered series with the realization
plot(ans, at.idx = NA, att.idx = 1, CI = NA)
lines(a, lty = "dotted")

## Check whether the residuals are Gaussian
plot(ans, type = "resid.qq")

## Check for linear serial dependence through 'acf'
plot(ans, type = "acf")


## <--------------------------------------------------------------------------->
## Example 2: Local level model for the Nile's annual flow.
## <--------------------------------------------------------------------------->
## Transition equation:
## alpha[t+1] = alpha[t] + eta[t], eta[t] ~ N(0, HHt)          
## Measurement equation:
## y[t] = alpha[t] + eps[t], eps[t] ~  N(0, GGt)

y <- Nile
y[c(3, 10)] <- NA  # NA values can be handled

## Set constant parameters:
dt <- ct <- matrix(0) 
Zt <- Tt <- matrix(1)
a0 <- y[1]            # Estimation of the first year flow 
P0 <- matrix(100)     # Variance of 'a0'

## Estimate parameters:
fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                   GGt = var(y, na.rm = TRUE) * .5),
                 fn = function(par, ...)
                 -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                 yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                 Zt = Zt, Tt = Tt, check.input = FALSE)

## Filter Nile data with estimated parameters:
fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]),
               GGt = matrix(fit.fkf$par[2]), yt = rbind(y))

## Compare with the stats' structural time series implementation:
fit.stats <- StructTS(y, type = "level")

fit.fkf$par
fit.stats$coef

## Plot the flow data together with fitted local levels:
plot(y, main = "Nile flow")
lines(fitted(fit.stats), col = "green")
lines(ts(fkf.obj$att[1, ], start = start(y), frequency = frequency(y)), col = "blue")
legend("top", c("Nile flow data", "Local level (StructTS)", "Local level (fkf)"),
       col = c("black", "green", "blue"), lty = 1)




#cleanEx()
#nameEx("plot.fkf")
### * plot.fkf

#flush(stderr()); flush(stdout())

### Name: plot.fkf
### Title: Plotting fkf Objects
### Aliases: plot.fkf
### Keywords: hplot

### ** Examples

## Local level model for the treering width data.
## Transition equation:
## alpha[t+1] = alpha[t] + eta[t], eta[t] ~ N(0, HHt)          
## Measurement equation:
## y[t] = alpha[t] + eps[t], eps[t] ~  N(0, GGt)

y <- treering
y[c(3, 10)] <- NA  # NA values can be handled

## Set constant parameters:
dt <- ct <- matrix(0) 
Zt <- Tt <- matrix(1)
a0 <- y[1]            # Estimation of the first width
P0 <- matrix(100)     # Variance of 'a0'

## Estimate parameters:
fit.fkf <- optim(c(HHt = var(y, na.rm = TRUE) * .5,
                   GGt = var(y, na.rm = TRUE) * .5),
                 fn = function(par, ...)
                 -fkf(HHt = matrix(par[1]), GGt = matrix(par[2]), ...)$logLik,
                 yt = rbind(y), a0 = a0, P0 = P0, dt = dt, ct = ct,
                 Zt = Zt, Tt = Tt, check.input = FALSE)

## Filter Nile data with estimated parameters:
fkf.obj <- fkf(a0, P0, dt, ct, Tt, Zt, HHt = matrix(fit.fkf$par[1]),
               GGt = matrix(fit.fkf$par[2]), yt = rbind(y))

## Plot the width together with fitted local levels:
plot(y, main = "Treering data")
lines(ts(fkf.obj$att[1, ], start = start(y), frequency = frequency(y)), col = "blue")
legend("top", c("Treering data", "Local level"), col = c("black", "blue"), lty = 1)

## Check the residuals for normality:
plot(fkf.obj, type = "resid.qq")

## Test for autocorrelation:
plot(fkf.obj, type = "acf", na.action = na.pass)



### * <FOOTER>
###
#cat("Time elapsed: ", proc.time() - get("ptime", pos = 'CheckExEnv'),"\n")
#grDevices::dev.off()
###
### Local variables: ***
### mode: outline-minor ***
### outline-regexp: "\\(> \\)?### [*]+" ***
### End: ***
#quit('no')
