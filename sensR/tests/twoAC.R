library(sensR)

## Testing border line cases:
twoAC(c(5, 0, 15))
twoAC(c(5, 0, 15), stat = "Wald")
twoAC(c(0, 5, 15))
twoAC(c(5, 15, 0))
twoAC(c(5, 0, 0))
twoAC(c(0, 0, 15))
twoAC(c(0, 5, 15), stat = "Wald")
twoAC(c(5, 15, 0), stat = "Wald")
twoAC(c(5, 0, 0), stat = "Wald")
twoAC(c(0, 0, 15), stat = "Wald")

## Testing twoAC functions:
(fm1 <- twoAC(c(2, 2, 6)))
confint(fm1, type = "Wald")
confint(fm1)
confint(fm1, level = .90)
confint(fm1, level = 1 - 1e-3)
confint(fm1, parm = "d.prime", type = "Wald")
confint(fm1, parm = "d.prime")
pr1 <- profile(fm1, range = c(-2, 3))
pr1 <- profile(fm1)
confint(pr1)

pr1 <- profile(fm1, alpha = 1e-5)
par(mfrow = c(2,2))
plot(pr1)
plot(pr1, Log = FALSE, relative = TRUE)
plot(pr1, Log = TRUE, relative = TRUE)
plot(pr1, Log = TRUE, relative = FALSE)

## Testing profile, confint and plot.profile in boundary cases:
## case 1:
(fit1 <- twoAC(c(5, 0, 0)))
logLik(fit1)
pr1 <- profile(fit1, range = c(-5, 2))
plot(pr1)
confint(pr1)

## case 2:
(fit2 <- twoAC(c(0, 5, 0)))
logLik(fit2)
pr2 <- try(profile(fit2, range = c(-5, 5)), silent = TRUE)
stopifnot(class(pr2) == "try-error")

## case 3
(fit3 <- twoAC(c(0, 0, 5)))
logLik(fit3)
pr3 <- profile(fit3, range = c(-5, 5))
plot(pr3)
confint(pr3)

## case 4:
(fit4 <- twoAC(c(5, 4, 0)))
logLik(fit4)
pr4 <- profile(fit4, range = c(-5, 5))
plot(pr4)
confint(pr4)

## case 5:
(fit5 <- twoAC(c(0, 4, 5)))
pr5 <- profile(fit5, range = c(-5, 5))
plot(pr5)
confint(pr5)

## case 6:
(fit6 <- twoAC(c(4, 0, 4)))
pr6 <- profile(fit6, range = c(-3, 3))
plot(pr6)
confint(pr6)

#################################
## testing clm2twoAC:

if(require(ordinal)) {
    response <- gl(3,1)
    fit.clm <- clm(response ~ 1, weights = c(2, 2, 6), link = "probit")
    fit.clm2 <- clm2(response ~ 1, weights = c(2, 2, 6), link = "probit")

    clm2twoAC(fit.clm)
    clm2twoAC(fit.clm2)
}


