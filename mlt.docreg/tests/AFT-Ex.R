
library("mlt")
library("eha")
set.seed(29)

## ************** Exponential - AFT *********************************

### true dgp
rY <- function(n, ...) rexp(n, ...)
pY <- function(x, ...) pexp(x, ...)
dY <- function(x, ...) dexp(x, ...)

gf <- gl(3, 1)
g <- rep(gf, 100)
y <- rY(length(g), rate = (1:nlevels(g))[g])
mydata <- data.frame(y = y, g = g)

boxplot(y ~ g, data = mydata)

Bb <- log_basis(numeric_var("y", support = range(y)), ui = "increasing",
                remove_intercept = TRUE)
Bx <- as.basis(~ g, data = mydata)
m <- ctm(Bb, shifting = Bx, todist = "MinExtrVal")

## Estimate coefficients
coef(opt <- mlt(m, data = mydata, fixed = c("log(y)" = 1)))

coef(aft <- survreg(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata,
                    dist = "exponential"))

coef(cox <- coxph(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata))

coef(phreg <- phreg(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata, 
                    dist = "weibull", shape = 1))

## Compare standard errors
## MLT
sqrt(diag(vcov(opt)))[c("g2", "g3")]
## Cox
sqrt(diag(vcov(cox)))
## phreg
sqrt(diag(phreg$var))[c("g2", "g3")]

## Compare log-Likelihoods
logLik(aft)
logLik(opt)


## Use a Weibull-AFT for estimation (Weibull shape parameter should be nu = 1)

## Estimate coefficients
(cf <- coef(opt2 <- mlt(m, data = mydata)))
cf[-1] / cf[1]

coef(aft2 <- survreg(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata,
                     dist = "weibull"))

## Compare Weibull shape paramters
1 / cf[1]
aft2$scale

## Compare log-Likelihoods
logLik(opt2)
logLik(aft2)

sqrt(diag(vcov(opt2)))[c("g2", "g3")]
sqrt(diag(vcov(aft2)))[c("g2", "g3")]


## *************** Right-censored

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE)), g = g)
coef(opt <- mlt(m, data = mydata, fixed = c("log(y)" = 1)))

## Estimate coefficients 
coef(aft <- survreg(y ~ g, data = mydata, dist = "expo"))
coef(cox <- coxph(y ~ g, data = mydata))
coef(phreg <- phreg(y ~ g, data = mydata, dist = "weibull", shape = 1))

## Compare standard errors
## MLT
sqrt(diag(vcov(opt)))[c("g2", "g3")]
## Cox
sqrt(diag(vcov(cox)))
## phreg
sqrt(diag(phreg$var))[c("g2", "g3")]



## ************** Left-censored

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE), 
                     type = "left"), g = g)

## Estimate coefficients
coef(opt <- mlt(m, data = mydata,  fixed = c("log(y)" = 1)))

coef(aft <- survreg(y ~ g, data = mydata, dist = "expo"))

## Compare standard errors
## MLT
sqrt(diag(vcov(opt)))[c("g2", "g3")]
## phreg
sqrt(diag(phreg$var))[c("g2", "g3")]


try(coef(cox <- coxph(y ~ g, data = mydata)))
try(coef(phreg <- phreg(y ~ g, data = mydata, dist = "weibull", shape = 1)))



## *************** Intervall-censored
mydata <- data.frame(y = Surv(y, y + 1, sample(0:3, length(y), replace = TRUE),
                     type = "interval"), g = g)

coef(opt<- mlt(m, data = mydata, fixed = c("log(y)" = 1)))
coef(aft <- survreg(y ~ g, data = mydata, dist = "expo"))

## Compare standard errors
## MLT
sqrt(diag(vcov(opt)))[c("g2", "g3")]
## phreg
sqrt(diag(phreg$var))[c("g2", "g3")]


try(coef(cox <- coxph(y ~ g, data = mydata)))
try(coef(phreg <- phreg(y ~ g, data = mydata, dist = "weibull", shape = 1)))




## ************** Weibull - AFT *********************************

set.seed(196)

### true dgp
rY <- function(n, ...) rweibull(n, ...)
pY <- function(x, ...) pweibull(x, ...)
dY <- function(x, ...) dweibull(x, ...)

gf <- gl(3, 1)
g <- rep(gf, 100)
y <- rY(length(g), scale = (1:nlevels(g))[g], shape = 3)
mydata <- data.frame(y = y, g = g)

boxplot(y ~ g, data = mydata)

Bb <- log_basis(numeric_var("y", support = range(y)), ui = "increasing", 
                remove_intercept = TRUE)
Bx <- as.basis(~ g, data = mydata)
m <- ctm(Bb, shifting = Bx, todist = "MinExtrVal")

## Estimate coefficients

## PH-scale
(cf <- coef(opt <- mlt(m, data = mydata)))

(coef_cox <- coef(cox <- coxph(Surv(y, rep(TRUE, nrow(mydata))) ~ g, 
                 	          data = mydata)))

(coef_phreg <- coef(phreg <- phreg(Surv(y, rep(TRUE, nrow(mydata))) ~ g,
                                   data = mydata, dist = "weibull")))

## AFT-scale
coef(aft <- survreg(Surv(y, rep(TRUE, nrow(mydata))) ~ g, data = mydata,
                    dist = "weibull"))

cf[-1] / cf[1]
coef_cox * aft$scale
coef_phreg[c("g2", "g3")] * aft$scale

## Compare shape parameters
1 / cf[1]
1 / exp(coef_phreg[c("log(shape)")])
aft$scale

## Compare standard errors
sqrt(diag(vcov(opt)))[c("g2", "g3")]
sqrt(diag(vcov(cox)))
sqrt(diag(phreg$var))[c("g2", "g3")]

## Compare log-Likelihoods
logLik(aft)
logLik(opt)


## ************************* Right-censored

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE)), g = g)

## Estimate coefficients
(cf <- coef(opt <- mlt(m, data = mydata)))
cf[-1] / cf[1]

coef(aft <- survreg(y ~ g, data = mydata, dist = "weibull"))

coef_cox <- coef(cox <- coxph(y ~ g, data = mydata))
coef_cox * aft$scale

coefs_phreg <- coef(phreg <- phreg(y ~ g, data = mydata, dist = "weibull"))
coefs_phreg[c("g2", "g3")] * aft$scale

## Compare standard errors
sqrt(diag(vcov(opt)))[c("g2", "g3")]
sqrt(diag(vcov(cox)))
sqrt(diag(phreg$var))[c("g2", "g3")]

## Compare log-Likelihoods
logLik(opt)
logLik(aft)


## ****************** Left-censored

mydata <- data.frame(y = Surv(y, sample(0:1, length(y), replace = TRUE), 
                     type = "left"), g = g)

## Estimate coefficients
(cf <- coef(opt <- mlt(m, data = mydata)))
cf[-1] / cf[1]

coef(aft <- survreg(y ~ g, data = mydata, dist = "weibull"))

## Compare log-Likelihoods
logLik(opt)
logLik(aft)


## ************** Interval-censored

mydata <- data.frame(y = Surv(y, y + 1, sample(0:3, length(y), replace = TRUE),
                     type = "interval"), g = g)

## Estimate coefficients
(cf <- coef(opt <- mlt(m, data = mydata)))
cf[-1] / cf[1]

coef(aft <- survreg(y ~ g, data = mydata, dist = "weibull"))

## Compare log-Likelihoods
logLik(opt)
logLik(aft)

