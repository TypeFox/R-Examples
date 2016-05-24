# Compound decisions for Gaussian Means and Variances  

n <- 100
r <- 10
m <- 2*r + 1
theta <- rep(c(1,2),n/2)
mu <- rep(c(2,3),each = n/2)
y <- rnorm(n*m, mean = rep(mu, each = m), sd = rep(theta, each = m))
id <- rep(1:n,each = m)
f <- WGLVmix(y,id, verb = 5)
require(lattice)
g <- expand.grid(theta = f$v, alpha = f$u)
g$fuv <- f$fuv
pl <- cloud(fuv ~ alpha * theta, data = g, type = "h", lwd = 2, 
      zlim = c(0, max(g$fuv)), scales = list(arrows = FALSE,
      xlab = expression(alpha), ylab = expression(theta), zlab = "density",
      screen = list(z = 10, x = -70)))
print(pl)
