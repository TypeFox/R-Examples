# syn.ppfm.r   Synthetic Poisson-Poisson Finite Mixture model
# Table 13.2: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
#
library(flexmix)
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
xb1 <- 1 + 0.25*x1 - 0.75*x2
xb2 <- 2 + 0.75*x1 - 1.25*x2
exb1 <- exp(xb2)
exb2 <- exp(xb1)
py1 <-  rpois(nobs, exb2)
py2 <-  rpois(nobs, exb1)
poixpoi <- py2
poixpoi <- ifelse(runif(nobs) > .9, py1, poixpoi)
pxp <- flexmix(poixpoi ~ x1 + x2, k=2, 
   model=FLXMRglm(family="poisson"))
summary(pxp)
parameters(pxp, component=1, model=1)
parameters(pxp, component=2, model=1)











library(MASS)     
library(gamlss)
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
xb <-  0.5 + 1.25*x1 - 1.5*x2     
delta <- .5      # value assigned to delta
exb <-exp(xb)
idelta <- (1/delta)*exb    
xg <-rgamma(50000, idelta, idelta, 1/idelta)
xbg <- exb*xg
nb1y <- rpois(50000, xbg)
jhgm <- gamlss(nbly ~ x1 + x2,family=NBII)
summary(jhgm)
# nb1 <- ml.nb1(los ~ hmo + white + type) # Table 10.8
# nb1







