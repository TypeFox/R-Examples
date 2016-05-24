## =============================================================================
## Set start values for parametric model fitting as a data frame.
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================


library("growthrates")
library("lattice")

## load data
data(bactgrowth)
splitted.data <- multisplit(value ~ time | strain + conc + replicate, data = bactgrowth)
dat <- splitted.data[[23]]

## initial parameters and box constraints
p   <- c(y0 = 0.03, mumax = .1, K = 0.1, alpha = 1, lambda = 2)
lower   <- c(y0 = 0.001, mumax = 1e-2, K = 0.005, alpha = -100, lambda = -20)
upper   <- c(y0 = 0.1,   mumax = 1,    K = 0.5,   alpha = 200,  lambda = 20)

## fit model
fit <- fit_growthmodel(FUN = grow_huang, p = p, time = dat$time, y = dat$value,
                       lower = lower, upper = upper,
                       control = list(trace = TRUE))

## coefficients and plot
coef(fit)
plot(fit)

## create data frame of parameters
ndata <- length(splitted.data)
pp <- as.data.frame(matrix(rep(p, ndata), byrow = TRUE, ncol = length(p)))
names(pp) <- names(p)

## fit 66 and 68 (T:62.5:2; R:125:2) need different start parameters
pp[66, ] <- pp[68, ] <- c(y0 = 0.01, mumax = .1, K = 0.02, alpha = 1, lambda = 2)


## fit growth models to all data
L <- all_growthmodels(value ~ time | strain + conc + replicate, FUN = grow_huang,
                      data = bactgrowth,
                      p = pp, lower = lower, upper = upper,
                      log = "y", ncores = 4, method = "L-BFGS-B")

par(mfrow = c(4,3))
plot(L)

res <- results(L)
