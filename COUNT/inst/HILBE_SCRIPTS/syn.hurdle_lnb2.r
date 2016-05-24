# syn.hurdle_lnb2.r  snthetic logit-NB2 hurdle model
# Table 11.4  Hilbe,JM (2011), Negative Binomial Regression, 2 ed, Cambridge Univ Press
library(MASS)
library(pscl)
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
xb <- 2 + .75*x1 - 1.25*x2
a <- .5
ia <- 1/.5
exb <- exp(xb)
xg <- rgamma(n = nobs, shape = a, scale = a)
xbg <-exb*xg
nby <- rpois(nobs, xbg)
nbdata <- data.frame(nby, x1, x2)
nby <- nbdata[nbdata$nby!=0, ]
pi <- 1/(1+exp(-(.9*x1 + .1*x2 + .2)))
bern <- runif(nobs)>pi
bern <- as.numeric(bern)
jhObs <- which( nbdata$bern==0 ) #  nbdata$nby <- ifelse(bern==0, 0, nbdata$nby)
nbdata$nby[jhObs] <- 0           #  hy <- nbdata$nby
hy <- nby
hlnb2 <- hurdle(hy ~ x1 + x2, dist="negbin", 
      zero.dist= "binomial", link="logit", data=nbdata)
summary(hlnb2)
