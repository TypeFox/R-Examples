## ----opts, echo = FALSE, message = FALSE---------------------------------
library("knitr")
#knitr::opts_chunk$set(eval = FALSE)

## ----eval=TRUE, echo=FALSE, results="hide"-------------------------------
suppressMessages(require("growthrates"))
#require("growthrates")

## ----eval=FALSE----------------------------------------------------------
#  library("growthrates")

## ----eval=TRUE-----------------------------------------------------------
data(bactgrowth)
str(bactgrowth)

## ----eval=TRUE-----------------------------------------------------------
head(bactgrowth)

## ---- fig.width=10, fig.height=14----------------------------------------
library(lattice)
data(bactgrowth)
xyplot(value ~ time|strain+as.factor(conc), data = bactgrowth, 
       groups = replicate, pch = 16, cex = 0.5)

## ------------------------------------------------------------------------
splitted.data <- multisplit(bactgrowth, c("strain", "conc", "replicate"))
dat <- splitted.data[[1]]

## ------------------------------------------------------------------------
fit <- fit_easylinear(dat$time, dat$value)

## ------------------------------------------------------------------------
summary(fit)

## ------------------------------------------------------------------------
coef(fit)      # exponential growth parameters
rsquared(fit)  # coefficient of determination (of log-transformed data)
deviance(fit)  # residual sum of squares of log-transformed data

## ---- fig.width=7--------------------------------------------------------
par(mfrow = c(1, 2))
plot(fit, log = "y")
plot(fit)

## ------------------------------------------------------------------------
fitx <- fit_easylinear(dat$time, dat$value, h = 8, quota = 0.95)
plot(fit)
lines(fitx, pch = "+", col = "blue")

## ---- fig.width=7--------------------------------------------------------
p     <- c(y0 = 0.01, mumax = 0.2, K = 0.1)
lower <- c(y0 = 1e-6, mumax = 0,   K = 0)
upper <- c(y0 = 0.05, mumax = 5,   K = 0.5)

fit1 <- fit_growthmodel(FUN = grow_logistic, p = p, dat$time, dat$value,
                        lower = lower, upper = upper)

p     <- c(yi = 0.02, ya = 0.001, kw = 0.1,	mumax = 0.2, K = 0.1)
lower <- c(yi = 1e-6, ya = 1e-6, kw = 0,    mumax = 0,   K = 0)
upper <- c(yi = 0.05, ya = 0.05, kw = 10,   mumax = 5,   K = 0.5)

fit2 <- fit_growthmodel(FUN = grow_twostep, p = p, time = dat$time, y = dat$value,
                        lower = lower, upper = upper)

coef(fit1)
coef(fit2)

par(mfrow = c(1, 2))
plot(fit1)
lines(fit2, col = "red")

plot(fit1, log = "y")
lines(fit2, col = "red")

## ------------------------------------------------------------------------
fit3 <- fit_growthmodel(FUN = grow_twostep, p = p, time = dat$time, y = dat$value,
                        lower = lower, upper = upper, which = c("kw", "mumax", "K"))

summary(fit3)

coef(fit3)
plot(fit3)

## ---- fig.width=7--------------------------------------------------------

dat <- splitted.data[[2]]
time <- dat$time
y    <- dat$value

## automatic smoothing with cv
res <- fit_spline(time, y)

par(mfrow = c(1, 2))
plot(res, log = "y")
plot(res)
coef(res)


## ----fig.width=14, fig.height=20-----------------------------------------
many_spline_fits <- all_splines(value ~ time | strain + conc + replicate,
                                data = bactgrowth, spar = 0.5)

par(mfrow = c(12, 6))
par(mar = c(2.5, 4, 2, 1))
plot(many_spline_fits)

## ------------------------------------------------------------------------
## initial parameters and box constraints
p   <- c(y0 = 0.03, mumax = .1, K = 0.1, h0 = 1)

lower   <- c(y0 = 0.001, mumax = 1e-2, K = 0.005, h0 = 0)
upper   <- c(y0 = 0.1,   mumax = 1,    K = 0.5,   h0 = 10)

## fit growth models to all data using log transformed residuals
many_baranyi1 <- all_growthmodels(
                   value ~ grow_baranyi(time, parms) | strain + conc + replicate,
                   data = bactgrowth,
                   p = p, lower = lower, upper = upper,
                   log = "y", ncores = 2)

## ------------------------------------------------------------------------
## use coefficients of first fit as new initial parameters
pp   <- coef(many_baranyi1)
## but set h0 to a fixed value
pp[, "h0"] <- 0.65
## re-fit models
many_baranyi2 <- all_growthmodels(
                   value ~ grow_baranyi(time, parms) | strain + conc + replicate,
                   data = bactgrowth,
                   p = pp, lower = lower, upper = upper,
                   which = c("y0", "mumax", "K"), log = "y", ncores = 2)

## ----fig.width=14, fig.height=20-----------------------------------------
par(mfrow = c(12, 6))
par(mar = c(2.5, 4, 2, 1))
plot(many_baranyi2)

## ---- fig.width=7, fig.height=3------------------------------------------
many_spline_res   <- results(many_spline_fits)
many_baranyi2_res <- results(many_baranyi2)
xyplot(mumax ~ log(conc+1)|strain, data = many_spline_res, layout = c(3, 1))
xyplot(mumax ~ log(conc+1)|strain, data = many_baranyi2_res, layout = c(3, 1))

