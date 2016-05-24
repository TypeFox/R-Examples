## =============================================================================
## This script demonstrates 'masking' of parameters,
## i.e. splitting the parameter set in fitted and fixed parameters
##
## Author: Thomas Petzoldt, TU Dresden
## License: GPL >= 2, https://www.gnu.org/licenses/
## Please cite our work when using this package.
## =============================================================================


library("growthrates")


data(bactgrowth)
splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))

## select a single data set
dat <- splitted.data[[7]]


p     <- c(yi=0.01, ya=0.01, kw=0.1,	mumax=0.2, K=0.1)
lower <- c(yi=1e-6, ya=1e-6, kw=0,    mumax=0,   K=0)
upper <- c(yi=0.05, ya=0.05, kw=10,   mumax=5,   K=0.5)


## fit of all parameters
fit1 <- fit_growthmodel(FUN=grow_twostep, p=p, time=dat$time, y=dat$value,
                        lower=lower, upper=upper,
                        method="L-BFGS-B")

p     <- c(yi=0.031, ya=0.00, kw=0.01,	mumax=0.2, K=0.1)
lower <- c(yi=1e-6, ya=1e-6, kw=0,    mumax=0,   K=0)
upper <- c(yi=0.05, ya=0.05, kw=10,   mumax=5,   K=0.5)

## fit only kw, mumax and K, but fix y_i and y_a
fit2 <- fit_growthmodel(FUN=grow_twostep, p=p, time=dat$time, y=dat$value,
                        lower=lower, upper=upper,
                        which = c("kw", "mumax", "K"),
                        method="L-BFGS-B")


p     <- c(y0=0.03, mumax=.176, K=.101, alpha=1, beta=1, gamma=1)
lower <- c(y0=0.01, mumax=0,   K=0, alpha=0.5, beta=0.5, gamma=1)
upper <- c(y0=0.04, mumax=5,   K=10, alpha=5, beta=10, gamma=10)


## generalized logistic
fit3 <- fit_growthmodel(FUN=grow_genlogistic, p=p, time=dat$time, y=dat$value,
                        lower=lower, upper=upper,
                        which = c("mumax", "K"),
                        method="L-BFGS-B")


p     <- c(y0=0.03, mumax=.176, K=.101, alpha=1, beta=1, gamma=1)
fit3 <- fit_growthmodel(FUN=grow_genlogistic, p=p, time=dat$time, y=dat$value,
                        lower=lower, upper=upper,
                        which = c("alpha", "beta", "gamma"),
                        method="L-BFGS-B")

coef(fit1)
coef(fit2)
coef(fit3)

plot(fit1)
lines(fit2, col="blue")
lines(fit3, col="red")

