
library(RcppZiggurat)

zsetseed(123456890)
print(zrnorm(10), digits=5)

zsetseed(123456890)
summary(zrnorm(1e6))
