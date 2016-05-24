# Likelihoods for Gaussian Mixtures of Means and Variances  

require(Rmosek)
require(gmler)

n <- 100
r <- 10
m <- 2*r + 1
theta <- rep(c(1,2),n/2)
mu <- rep(c(2,3),each = n/2)
y <- rnorm(n*m, mean = rep(mu, each = m), sd = rep(theta, each = m))
id <- rep(1:n,each = m)
f <- GLVmix(y,id)
X11(width = 8, height = 5)
par(mfrow = c(1,2))
plot(f$u,f$fu,type = "l")
plot(f$v,f$fv,type = "l")
ha <- f$logLik
X11(width = 8, height = 9)
par(mfrow = c(5,5))
vs <- 16:40/10
h <- vs
for(i in 1:length(vs)){
	g <- GLVmix(y,id, v = vs[i], verb = 5)
	plot(g$u, g$fu, type = "l", xlab = "", main = paste("s = ", vs[i]))
	h[i] <- g$logLik
	}
X11(width = 6, height = 6)
plot(vs, h, xlab = expression(sigma^2), ylab = "log likelihood")
lines(vs,h)
