# syn.geo.r     Synthetic log-geometric regression
# Table 10.4: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
library(MASS)
nobs <- 50000
x2 <- runif(nobs)
x1 <- runif(nobs)
xb <- 2*x1 - .5*x2 - 1   
exb <- exp(xb)               
xg <- rgamma(n = nobs, shape = 1, rate = 1) 
xbg <-exb*xg                 
gy <- rpois(nobs, xbg)      
gnb2 <-glm.nb(gy ~ x1 + x2)   
summary(gnb2)
gpy <- glm(gy ~ x1 + x2, family=poisson)
summary(gpy)




