# Demo for the maximum likelihood estimation of parameters from
# some selected distributions
# At the moment this is copied from some .Rd file

## Negative binomial distribution
## Data from Bliss and Fisher (1953).

appletree <- data.frame(y = 0:7, w = c(70, 38, 17, 10, 9, 3, 2, 1))
fit <- vglm(y ~ 1, negbinomial(deviance = TRUE), data = appletree,
            weights = w, crit = "coef", half.step = FALSE)
summary(fit)
coef(fit, matrix = TRUE)
Coef(fit)
deviance(fit)  # NB2 only; needs 'crit = "coef"' & 'deviance = TRUE' above


## Beta distribution

set.seed(123)
bdata <- data.frame(y = rbeta(nn <- 1000, shape1 = exp(0), shape2 = exp(1)))
fit1 <- vglm(y ~ 1, betaff, data = bdata, trace = TRUE)
coef(fit1, matrix = TRUE)
Coef(fit1)  # Useful for intercept-only models

# General A and B, and with a covariate
bdata <- transform(bdata, x2 = runif(nn))
bdata <- transform(bdata, mu   = logit(0.5 - x2, inverse = TRUE),
                          prec =   exp(3.0 + x2))  # prec == phi
bdata <- transform(bdata, shape2 = prec * (1 - mu),
                          shape1 = mu * prec)
bdata <- transform(bdata,
                   y = rbeta(nn, shape1 = shape1, shape2 = shape2))
bdata <- transform(bdata, Y = 5 + 8 * y)  # From 5 to 13, not 0 to 1
fit2 <- vglm(Y ~ x2, data = bdata, trace = TRUE,
             betaff(A = 5, B = 13, lmu = elogit(min = 5, max = 13)))
coef(fit2, matrix = TRUE)




