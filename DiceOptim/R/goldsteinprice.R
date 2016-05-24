goldsteinprice <- function(x)
{
# Goldstein and Price test function (standardized version)
# --------------------------------------------------------
# Dimension: n = 2
# Number of local minima: 4
# The global minimum: 
# x* = c(0.5, 0.25), f(x*) = -3.129172
# The local minima:
# x*,2 = c(0.35, 0.4),  f(x*,2) = -2.180396
# x*,3 = c(0.95, 0.55), f(x*,3) = -1.756143
# x*,4 = c(0.8 , 0.7),  f(x*,4) = -0.807367

m <- 8.6928
s <- 2.4269

x1 <- 4 * x[1] - 2
x2 <- 4 * x[2] - 2

a <- 1 + (x1+x2+1)^2 * (19 - 14*x1 + 3*x1^2 - 14*x2 + 6*x1*x2 + 3*x2^2)
b <- 30 + (2*x1 - 3*x2)^2 * (18 - 32*x1 + 12*x1^2 + 48*x2 - 36*x1*x2 + 27*x2^2)
f <- log(a*b)
f <- (f-m)/s

return(f)
}