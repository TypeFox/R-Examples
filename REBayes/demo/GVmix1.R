# Compound decisions for Gaussian Variances  

# Initialize Rmosek
require(Rmosek)

n <- 100
r <- 10
m <- 2*r + 1
theta <- rep(c(.2,.3),n/2)
y <- rnorm(n*m, mean = 0, sd = rep(theta, each = m))
id <- rep(1:n,each = m)
x <- tapply(y, id, var)
m <- tapply(y, id, length)
f <- GVmix(x,m,  verb = 0) 
g <- GVmix(x,m, rtol = 1e-10, verb = 0) # Needs stricter rtol
plot(g,xlab = expression(sigma^2), main = "Estimated Mixing Density")
lines(g, col = 4)
abline(v = .04, col = "red")
abline(v = .09, col = "red")
