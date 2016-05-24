# nbr2_7_1.r   
# Table 7.1: Hilbe, Negative Binomial Regression, 2 ed, Cambridge Univ Press 
# Poisson with frequency table of observed counts

library(MASS)         
nobs <- 50000
x1 <- runif(nobs)
x2 <- runif(nobs)
x3 <- runif(nobs)
py <-rpois(nobs, exp(1 + 0.5*x1 - 0.75*x2 + 0.25*x3))
cnt <-  table(py)
df <- data.frame( prop.table( table(py) ) )
df$cumulative <- cumsum( df$Freq )
dfall <- data.frame(cnt, df$Freq, df$cumulative)
dfall







