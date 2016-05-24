## =============================================================================
## Determine growth rates with parametric models given in closed form or as
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
names(splitted.data)
dat <- splitted.data[["D:0.49:1"]]

## Logistic growth model
p     <- c(y0 = 0.01, mumax = 0.03, K = 0.1)
lower <- c(y0 = 1e-6, mumax = 0,   K = 0)
upper <- c(y0 = 0.05, mumax = 5,   K = 0.5)

fit1 <- fit_growthmodel(FUN = grow_logistic, dat$time, dat$value, p = p,
                        control = list(trace = TRUE),
                        lower = lower, upper = upper)

fit2 <- fit_growthmodel(FUN = grow_logistic, dat$time, dat$value, p = p,
                        lower = lower, upper = upper,
                        transform = "log", control = list(trace = TRUE))


## Two step model with lag phase given as system of differential equations
p     <- c(yi = 0.01, ya = 0.01, kw = 0.1,	mumax = 0.2, K = 0.1)
lower <- c(yi = 1e-6, ya = 1e-6, kw = 0,    mumax = 0,   K = 0)
upper <- c(yi = 0.05, ya = 0.05, kw = 10,   mumax = 5,   K = 0.5)


fit3 <- fit_growthmodel(FUN = grow_twostep, p = p, time = dat$time, y = dat$value,
                        lower = lower, upper = upper,
                        control = list(trace = TRUE),
                        method = "L-BFGS-B")

fit4 <- fit_growthmodel(FUN = grow_twostep, dat$time, dat$value, p = p,
                        lower = lower, upper = upper, transform = "log",
                        control = list(trace = TRUE),
                        method = "L-BFGS-B")

## show results graphically
plot(fit1)#, log = "y")
lines(fit2, col = "red")
lines(fit3, col = "blue", lty = "dotted")
lines(fit4, col = "cyan", lty = "dashed")


## Exercise:
##   try to improve the model formulation to get a better fit.
##   see example user_defined_ode_model.R

