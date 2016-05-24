# syn.poissono.r   Poisson with offset
# Table 6.20 : Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
off <- rep(1:5, each=10000, times=1)*100 # offset 
loff <- log(off)                         
py <-rpois(nobs, exp(2 + .75*x1 -1.25*x2 + loff))  
poir <-glm(py ~ x1 + x2 + offset(loff), family=poisson)
summary(poir)
confint(poir)







