# syn.nb1.r     Synthetic NB1 regression
# Table 10.6: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
# Synthetic NB1 regression   amend as needed
#
library(MASS)     
library(gamlss)
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
xb <-  0.5 + 1.25*x1 - 1.5*x2     
delta <- .5      # value assigned to delta
exb <-exp(xb)
idelta <- (1/delta)*exb    
xg <-rgamma(n = 50000, shape = idelta, rate = idelta)
xbg <- exb*xg
nb1y <- rpois(50000, xbg)
jhgm <- gamlss(nbly ~ x1 + x2,family=NBII)
summary(jhgm)
# nb1 <- ml.nb1(los ~ hmo + white + type) # Table 10.8
# nb1







