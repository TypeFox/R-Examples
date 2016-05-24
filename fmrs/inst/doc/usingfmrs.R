## ------------------------------------------------------------------------
library(fmrs)
set.seed(1980)
nComp = 2
nCov = 10
nObs = 500
dispersion = c(1, 1)
mixProp = c(0.4, 0.6)
rho = 0.5
coeff1 = c( 2,  2, -1, -2, 1, 2, 0, 0,  0, 0,  0)
coeff2 = c(-1, -1,  1,  2, 0, 0, 0, 0, -1, 2, -2)
umax = 40

## ------------------------------------------------------------------------
dat <- fmrs.gendata(nObs = nObs, nComp = nComp, nCov = nCov,
                     coeff = c(coeff1, coeff2), dispersion = dispersion,
                     mixProp = mixProp, rho = rho, umax = umax, 
                     disFamily = "lnorm")

## ------------------------------------------------------------------------
res.mle <- fmrs.mle(y = dat$y, x = dat$x, delta = dat$delta,
                    nComp = nComp, disFamily = "lnorm",
                    initCoeff = rnorm(nComp*nCov+nComp),
                    initDispersion = rep(1, nComp),
                    initmixProp = rep(1/nComp, nComp))

coefficients(res.mle)
dispersion(res.mle)
mixProp(res.mle)

## ------------------------------------------------------------------------
res.lam <- fmrs.tunsel(y = dat$y, x = dat$x, delta = dat$delta,
                       nComp = nComp, disFamily = "lnorm",
                       initCoeff = coefficients(res.mle),
                       initDispersion = dispersion(res.mle),
                       initmixProp = mixProp(res.mle), 
                       penFamily = "adplasso")
show(res.lam)

## ------------------------------------------------------------------------
res.var <- fmrs.varsel(y = dat$y, x = dat$x, delta = dat$delta,
                       nComp = nComp, disFamily = "lnorm",
                       initCoeff = coefficients(res.mle),
                       initDispersion = dispersion(res.mle),
                       initmixProp = mixProp(res.mle), 
                       penFamily = "adplasso",
                       lambPen = slot(res.lam, "lambPen"))
coefficients(res.var)
dispersion(res.var)
mixProp(res.var)

## ------------------------------------------------------------------------
round(coefficients(res.var)[-1,],5)

