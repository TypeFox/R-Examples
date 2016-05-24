## =============================================================================
## Determine growth rates with a heuristic linear method, similar to
##   the method of Hall et al. 2013,  doi:10.1093/molbev/mst197
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================


library("growthrates")
library(lattice)

data(bactgrowth)

## select a single data set
dat <- splitted.data[[1]]
plot(value ~ time, data = dat)

## fit linear model to steepest part of the log-transformed data
fit <- fit_easylinear(dat$time, dat$value)
plot(fit)
plot(fit, log = "y")

## diagnostic plot
plot(fit, which = "diagnostics")

## change settings of the algorithm
fitx <- fit_easylinear(dat$time, dat$value, h = 8, quota = 0.95)

par(mfrow=c(1, 2))
plot(fit, log="y")
lines(fitx, pch = "+", col = "blue")

plot(fit)
lines(fitx, pch = "+", col = "blue")

## fit all data sets
allFits <- all_easylinear(value ~ time|strain + conc + replicate,
                          data = bactgrowth)
coef(allFits)
