cat("\ntest corMatern:")

data(blackcap)
blackcapD <-cbind(blackcap,dummy=1) ## obscure, isn't it? 
require(nlme)
## With method= 'ML' in lme, The correlated random effect is described 
##  as a correlated residual error and no extra residual variance is fitted:
bf <- lme(fixed = migStatus ~ means, data = blackcapD, random = ~ 1 | dummy, 
    correlation = corMatern(form = ~ longitude+latitude  | dummy), 
    method = "ML")

expect_equal(logLik(bf)[[1]],-7.941674,tolerance=1e-6)
expect_equal(exp((bf$modelStruct$corStruct)[[1]]),18.35958,tolerance=1e-5) ## nu
expect_equal(exp((bf$modelStruct$corStruct)[[2]]),0.6285744,tolerance=1e-5) ## range
