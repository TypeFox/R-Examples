if(!require("GNE"))stop("this test requires package GNE.")

GrAphiFB(0, 1)
GrAphiFB(0, 0)
library(GNE)
GrAphiFB(0, 0)
GrAphiFB(0, 1)
GrAphiFB(0, 1/10)
GrAphiFB(1/10, 0)
GrAphiFB(0:10, 0:10)
0:10 == 0 & 1:10 == 0
0:10 == 0 & 0:10 == 0
a <- 0:10
b <- 0:10
a / sqrt(a^2+b^2) - 1
a / sqrt(a^2+b^2) - 1
a <- c(0, rnorm(10))
b <- c(0, rnorm(10))
a
b
a / sqrt(a^2+b^2) - 1
GrAphiFB(a, b)
GrBphiFB(a, b)
a <- cbind(c(0, rnorm(10)), rnorm(11))
a
b <- cbind(c(0, rnorm(10)), rnorm(11))
GrAphiFB(a, b)
a == 0 & b == 0
a / sqrt(a^2+b^2) - 1
GrBphiFB(a, b)


xmax <- 10
x <- seq(-xmax, xmax, length=31)
p <- 0

sapply(x, function(x) phipFB(x, 1, p) - phiFB(x, 1))

p <- pi

sapply(x, function(x) phipFB(x, 1, p) - (sqrt(x^2+1) - (x+1) - p*pmax(x, 0)*max(1,0)))
sapply(x, function(x) phipFB(x, -1, p) - (sqrt(x^2+1) - (x-1) - p*pmax(x, 0)*max(-1,0)))

