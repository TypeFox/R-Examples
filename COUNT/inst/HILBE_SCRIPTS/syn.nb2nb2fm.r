# syn.nb2nb2fm.r   Synthetic NB2-NB2 Finite Mixture model
# Table 13.3: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
library(gamlss.mx)
nobs <- 50000
x1 <- (runif(nobs))      
x2 <- qnorm(runif(nobs))      
xb1 <- 1 + .25*x1 - .75*x2
xb2 <- 2 + .75*x1 - 1.25*x2     
a1 <- .5                        
a2 <- 1.5
ia1 <- 1/a1
ia2 <- 1/a2                     
exb1 <- exp(xb1) 
exb2 <- exp(xb2)                 
xg1 <- rgamma(n = nobs, shape = a1, rate = a1)  
xg2 <- rgamma(n = nobs, shape = a2, rate = a2)  
xbg1 <-exb1*xg1 
xbg2 <-exb2*xg2                  
nby1 <- rpois(nobs, xbg1)
nby2 <- rpois(nnobs, xbg2)       
nbxnb <- nby2
nbxnb <- ifelse(runif(nobs) > .9, nby1, nbxnb)
nxn <- gamlssNP(nbxnb~x1+x2, random=~1,family=NBI, K=2)
summary(nxn)

