library(tsDyn)

mod.lstar <- lstar(log10(lynx), m=2, mTh=c(0,1), control=list(maxit=3000))
mod.lstar
deviance(mod.lstar)
c(AIC(mod.lstar),BIC(mod.lstar))

mod.lstar2 <- lstar(log10(lynx), m=1, control=list(maxit=3000))
mod.lstar2
deviance(mod.lstar2)
c(AIC(mod.lstar2),BIC(mod.lstar2))

## include: none
mod.lstar_noConst <- lstar(log10(lynx), m=2, control=list(maxit=1000), include="none")
mod.lstar_noConst
deviance(mod.lstar_noConst)
c(AIC(mod.lstar_noConst),BIC(mod.lstar_noConst))

## include: trend
mod.lstar_trend <- lstar(log10(lynx), m=2, control=list(maxit=1000), include="trend")
mod.lstar_trend
deviance(mod.lstar_trend)
c(AIC(mod.lstar_trend),BIC(mod.lstar_trend))

## include: both
mod.lstar_both <- lstar(log10(lynx), m=2, control=list(maxit=1000), include="both")
mod.lstar_both
deviance(mod.lstar_both)
c(AIC(mod.lstar_both),BIC(mod.lstar_both))

## grid attributes
mod.lstar3 <- lstar(log10(lynx), m=2, control=list(maxit=3000), starting.control=list(gammaInt=c(1,1000), nTh=100))
mod.lstar3
deviance(mod.lstar3)
c(AIC(mod.lstar3),BIC(mod.lstar3))


mod.lstar_ALL <- list(mod.lstar=mod.lstar, mod.lstar2=mod.lstar2, 
                      mod.lstar_noConst=mod.lstar_noConst,mod.lstar_trend=mod.lstar_trend,
                      mod.lstar_both=mod.lstar_both,mod.lstar3=mod.lstar3)

sapply(mod.lstar_ALL, function(x) c(AIC=AIC(x), BIC=BIC(x), deviance=deviance(x)))
sapply(mod.lstar_ALL, function(x) tail(coef(x),4))
sapply(mod.lstar_ALL, function(x) tail(coef(x,hyperCo=FALSE),4))
sapply(mod.lstar_ALL, function(x) head(x$model,2))
