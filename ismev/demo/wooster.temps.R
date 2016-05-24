data(wooster)
x <- seq(along = wooster)
usin <- function(x, a, b, d)
  a + b * sin(((x-d) * 2 * pi)/365.25)
wu <- usin(x, -30, 25, -75)
ydat <- cbind(sin(2*pi*x/365.25), cos(2*pi*x/365.25))
wooster.pp <- pp.fit(-wooster, threshold = wu, ydat = ydat, mul = 1:2, sigl = 1:2, siglink = exp, method = "BFGS")
pp.diag(wooster.pp)

