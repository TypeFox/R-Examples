# syn.nb2o.r     Synthetic NB2 with offset
# Table 9.9: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
library(MASS)
x1 <- qnorm(runif(50000))
x2 <- qnorm(runif(50000))
off <- rep(1:5, each=10000, times=1)*100  # offset 
loff <- log(off)                          # log of offset
xb<-2 + .75*x1 -1.25*x2 + loff            # linear predictor
exb <-exp(xb)                  # inverse link
a <- .5                        # assign value to alpha
ia <- 1/.5                     # invert alpha
xg <- rgamma(n = 50000, shape = a, rate = a)  # generate gamma variates w alpha
xbg <-exb*xg                   # mix Poisson and gamma variates
nbyo <- rpois(50000, xbg)      # generate NB2 variates - w offset
nb2o <-glm.nb(nbyo ~ x1 + x2 + offset(loff))  # model NB2
summary(nb2o)








