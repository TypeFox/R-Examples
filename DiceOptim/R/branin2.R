branin2 <- function(x)
{
# Branin test function (standardized version)
#--------------------------------------------
# Dimension: n = 2
# Number of local minima: 3 (the global ones)
# The global minima: 
# x*,1 = c(0.1239, 0.8183), f(x*,1) = -1.047410
# x*,2 = c(0.5428, 0.1517), f(x*,2) = -1.047410
# x*.3 = c(0.9617, 0.1650), f(x*,3) = -1.047410

m <- 54.8104
s <- 51.9496

xx <- 15 * x[1] - 5
y <- 15 * x[2]

f <- (y - 5.1*xx^2/(4*pi^2) + 5*xx/pi - 6)^2 + 10*(1 - 1/(8*pi))*cos(xx) + 10
f <- (f-m)/s

return(f)
}
