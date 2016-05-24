# Demo for WTLVmix
n <- 100
r <- 10
m <- 2* r + 1
theta <- rep(c(1,4), n/2)
mu <- rep(c(2,3), each = n/2)
y <- rnorm(n * m, mean = rep(mu, each = m), sd = rep(sqrt(theta), each = m))
id <- rep(1:n, each = m)
f <- WTLVmix(y, id, verb = 5)
X11(width = 8, height = 5)
par(mfrow = c(1,2))
plot(f$u, f$fu, main = expression(paste("Density of ", alpha, sep = "")),
     xlab = expression(alpha), ylab = expression(f(alpha)), type="l")
plot(f$v, f$fv, main = expression(paste("Density of ", theta, sep = "")),
     xlab = expression(theta), ylab = expression(f(theta)), type="l")
