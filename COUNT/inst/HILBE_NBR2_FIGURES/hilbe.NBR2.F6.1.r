# hilbe.NBR2.F6.1.r
# Table 6.4 plus added code
# Synthetic Poisson model with a user defined values; 
#   graphic of predicted mean values
# From Hilbe, Negative Binomial regression, 2nd ed, Cambridge Univ. Press
# Table 6.4; Figure 6.1
#
nobs <- 50000
x1 <- qnorm(runif(nobs))
x2 <- qnorm(runif(nobs))
py <-rpois(nobs, exp(2 + .75*x1 -1.25*x2))
mpy <- mean(py)    
ypoi <- (exp(-mpy)*mpy^py)/gamma(py+1)
 plot(ypoi~ py, xlim=c(0,50), main="Synthetic Poisson Model: Mean=21")



