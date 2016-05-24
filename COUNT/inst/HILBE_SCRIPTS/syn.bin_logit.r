# syn.bin_logit.r     Synthetic grouped logistic regression
# Table 9.22: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
nobs <- 50000                                  
x1 <- runif(nobs)                      
x2 <- runif(nobs)                      
d <- rep(1:5, each=10000, times=1)*100 # denominator
xb <- 2 + .75*x1 - 1.25*x2             # linear predictor; values
exb <- 1/(1+exp(-xb))                                       # fit; predicted prob
by  <- rbinom(nobs, size = d, p = exb) # random binomial variate
dby <- d - by                          # denominator - numerator
gby <- glm(cbind(by,dby) ~ x1 + x2, family=binomial(link="logit"))
summary(gby)                           # displays model output










