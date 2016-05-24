# geo_rng.r     
# Table 10.2: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
# Synthetic geometric regression
#
library(MASS)
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
xb <- 2 + .75*x1 - 1.25*x2   # parameter values
exb <- exp(xb)               # Poisson predicted value
xg <- rgamma(n = nobs, shape = 1, rate = 1)  # gamma variate, param 1,1
xbg <-exb*xg                 # mix Poisson and gamma variates
gy <- rpois(nobs, xbg)       # generate NB2 variates
geo <-glm.nb(gy ~ x1 + x2)   # model geometric
summary(geo)
nbg <- glm(gy ~ x1 + x2,, family=negative.binomial(1))
summary(nbg)                 # GLM NB2 with ?=1    


