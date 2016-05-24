## =============================================================================
## Determine growth rates with a nonparametric smoothing spline.
##   system of differential equations.
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================

library("growthrates")

data(bactgrowth)
splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))

## select a single data set
dat <- splitted.data[[2]]
time <- dat$time
y    <- dat$value

## automatic smoothing with cv
res <- fit_spline(time, y)

plot(res, log="y")
plot(res)
coef(res)

## a more difficult data set
dat <- splitted.data[[56]]
time <- dat$time
y <- dat$value

## default parameters
res <- fit_spline(time, y)
plot(res, log="y")

## small optgrid, trapped in local minimum
res <- fit_spline(time, y, optgrid=5)
plot(res, log="y")


## manually selected smoothing parameter
res <- fit_spline(time, y, spar=.5)
plot(res, log="y")
plot(res, ylim=c(0.005, 0.03))

## todo
summary(res)
residuals(res) #  log-transformed!
rsquared(res)  #  log-transformed
deviance(res)  #

