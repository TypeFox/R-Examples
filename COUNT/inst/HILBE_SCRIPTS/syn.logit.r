# syn.logit.r     Synthetic logistic regression
# Table 9.20: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
library(MASS)
nobs <- 50000                           
x1 <- runif(nobs)                       
x2 <- runif(nobs)                       
xb <- 2 + .75*x1 - 1.25*x2                 # linear predictor 
exb <- 1/(1+exp(-xb))                      # fit; predicted prob
by  <- rbinom(nobs, size = 1, prob =exb)   # random logit variates
lry <- glm(by ~ x1 + x2, family=binomial(link="logit"))
summary(lry)                            









