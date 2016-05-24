# syn.poisson.r  
# Table 6.4 : Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
nobs <- 50000
x1 <- qnorm(runif(nobs))
x2 <- qnorm(runif(nobs))
py <-rpois(nobs, exp(2 + .75*x1 -1.25*x2))
jhpoi <-glm(py ~ x1 + x2, family=poisson)
summary(jhpoi)
confint(jhpoi)

