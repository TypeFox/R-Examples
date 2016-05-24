rosenbrock4 <- function(x)
{
# 4D-Rosenbrock test function (standardized version)
#---------------------------------------------------
# Dimension: n = 4
# Number of local minima: 2
# The global minimum: 
# x* = c(0.4,0.4,0.4,0.4), f(x*) = -1.019701
# The local minimum:
# x*,2 = c(0.26667,0.4,0.4,0.4), f(x*,2) = -1.019691

m <- 382658.057227524
s <- 375264.858362295

x <- 15 * x - 5

x1 <- x[seq(1,3)]
x2 <- x[seq(2,4)]
f <- sum(100*(x2 - x1^2)^2 + (1 - x1)^2);
f <- (f-m)/s
return(f)
}
