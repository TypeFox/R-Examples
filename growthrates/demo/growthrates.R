library("growthrates")

## =============================================================================
## single fits
## =============================================================================

data(bactgrowth)
splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))

## make subset of single experiment
dat <- splitted.data[["D:0:1"]]

fit1 <- fit_spline(dat$time, dat$value)
par(mfrow=c(1, 2))
plot(fit1, log="y")
plot(fit1)

## derive start parameters from spline fit
p <- coef(fit1)

## fit exponential to subset of first 10 data of the experiment
first10 <-  dat[1:10, ]
fit2    <- fit_growthmodel(grow_exponential, p=p, time=first10$time, y=first10$value)

## fit logistic curve to all data of the subset
p <- c(coef(fit1), K = max(dat$value))
fit3 <- fit_growthmodel(grow_logistic, p=p, time=dat$time, y=dat$value, transform="log")

plot(fit1)
lines(fit2, col="green")
lines(fit3, col="red")

## =============================================================================
## multiple fits
## =============================================================================

## use spline method with user-definded smoothness (spar)
sfit <- all_splines(value ~ time | strain + conc + replicate, data = bactgrowth,
                    spar=0.5)

## initial parameters
(p <- c(coef(fit1), K = max(dat$value)))

## avoid negative parameters
lower = c(y0=0, mu=0, K=0)

## total fit to all data
pfit1 <- all_growthmodels(value ~ grow_logistic(time, parms), data=bactgrowth,
                         p = p, lower = lower)

par(mfrow=c(1,1))
plot(pfit1)
## fit all models
pfit2 <- all_growthmodels(value ~ grow_logistic(time, parms) |
                           strain + conc + replicate, data=bactgrowth,
                         p = p, lower = lower, ncores=2)
par(mfrow=c(4, 3))
plot(pfit2)
