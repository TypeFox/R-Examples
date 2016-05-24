
set.seed(100)
n <- 100
p <- 20
x <- matrix(rnorm(n*p), n,p)
y <- x%*%c(rep(0.2,5),rep(0,p-5)) + rnorm(n)
eta <- x%*%c(rep(0.2,5),rep(0,p-5))
ppi <- 1/(1+exp(-eta))
yy <- rbinom(n, 1, ppi)
rbinom(n,1,0.5)
y <- yy

out <- cvplogistic(y, x, "mcp", "mmcd", 1/2.7)
cvplogistic(y, x, "mcp", "adaptive")
cvplogistic(y, x, "mcp", "llacda")
out1 <- cvplogistic(y, x, "scad", "mmcd", 0)

out2 <- cvplogistic(y, x, kappa=0)

hybrid.logistic(y, x, "mcp", 0)
hybrid.logistic(y, x, "scad", 0)

cv.cvplogistic(y, x, "mcp", "mmcd")
cv.cvplogistic(y, x, "mcp", "adaptive")
cv.cvplogistic(y, x, "mcp", "llacda")
cv.cvplogistic(y, x, "scad", "mmcd")

cv.hybrid(y,x, "mcp")
cv.hybrid(y,x, "scad")

path.plot(out)
