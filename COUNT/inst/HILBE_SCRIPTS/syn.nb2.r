# syn.nb2.r   Synthetic NB2
# Table 9.3: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
library(MASS)                 
nobs <- 50000
x1 <- qnorm(runif(nobs))      # random normal N[0,1] variate
x2 <- qnorm(runif(nobs))      # random normal N[0,1] variate
xb <- 2 + .75*x1 - 1.25*x2    # parameter values
a <- .5                       # assign value to ancillary parameter
ia <- 1/.5                    # invert alpha
exb <- exp(xb)                # Poisson predicted value
xg <- rgamma(n = nobs, shaep = a, rate = a)  # generate gamma variates given alpha
xbg <-exb*xg                  # mix Poisson and gamma variates
nby <- rpois(nobs, xbg)       # generate NB2 variates
jhnb2 <-glm.nb(nby ~ x1 + x2) # model NB2
summary(jhnb2)



