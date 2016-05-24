# syn.probit.r     Synthetic probit regression
# Table 9.21: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
library(MASS)
nobs <- 50000                           
x1 <- runif(nobs)                       
x2 <- runif(nobs)                       
xb <- 2 + .75*x1 - 1.25*x2                  
exb <- pnorm(xb)
by  <- rbinom(nobs, size = 1, prob =exb)
pry <- glm(by ~ x1 + x2, family=binomial(link="probit"))
summary(pry)








