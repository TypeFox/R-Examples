# syn.cgeo.r     Synthetic canonical geometric regression
# Table 10.3: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
library(MASS)
nobs <- 50000
x2 <- runif(nobs)
x1 <- runif(nobs)
xb <- 1.25*x1 + .1*x2 - 1.5   
mu <- 1/(exp(-xb)-1)
p <- 1/(1+mu)
r <- 1
gcy <- rnbinom(nobs, size=r, prob = p)
source("c://source/ml.nbc.r")  # NB-C function
g2y <- ml.nbc(y ~ x1 + x2)
summary(g2y)
library(gamlss)
hnbii <- gamlss(y~ x1 + x2, family=NBII)
summary(hnbii) 


