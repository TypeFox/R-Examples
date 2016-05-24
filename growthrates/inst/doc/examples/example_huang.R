## =============================================================================
## Test of the Huang Growth model, cf. Huang, L (2011)
## Int. J. Food Microbiology,  10.1016/j.fm.2010.05.019
##
## Note: the original model formulation works in log space, so that
##       y(t), y_0 and y_max are all given as natural log.
##  This is advantageous for model fitting, but here we use y_0 and K in
##  untransformed space to be compatible with the growth parameters of
##  most other models. The downside is, that we need box constraints in most
##  cases and possibly more iterations.
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================


library("growthrates")
library("lattice")

## load data
data(bactgrowth)
splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
dat <- splitted.data[[23]]

## initial parameters and bocx constraints
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


## fit growth models to all data using (log transformed residuals)
L <- all_growthmodels(value ~ grow_huang(time, parms) | strain + conc + replicate,
                      data = bactgrowth,
                      p = p, lower = lower, upper = upper,
                      log = "y")

par(mfrow = c(4,3))
plot(L, log = "y")

par(mfrow = c(4,3))
plot(L)

res <- results(L)
xyplot(lambda ~ log(conc + 1)| strain, data = res)
xyplot(alpha ~ log(conc + 1)| strain, data = res)

## and most importantly, the max growth rates
xyplot(mumax ~ log(conc + 1)| strain, data = res)

## 2nd approach: fit selected parameters, fix the remaining (here: alpha)
## alpha = 4 from the IPMP tutorial
p   <- c(y0 = 0.03, mumax = .1, K = 0.1, alpha = 4, lambda = 2)
L2 <- all_growthmodels(value ~ grow_huang(time, parms) | strain + conc + replicate,
                      data = bactgrowth,
                      p = p, lower = lower, upper = upper,
                      which = c("y0", "mumax", "K", "lambda"),
                      method = "Marq", log = "y")

par(mfrow = c(4,3))
plot(L2)

res2 <- results(L2)
xyplot(mumax ~ log(conc + 1)| strain, data = res2)

plot(res$mumax, res2$mumax)
