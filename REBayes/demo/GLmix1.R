# Ref.  Jiang and Zhang (Annals, 2009): 

# Test problem:  y ~ N(mu_i , 1), mu_i ~ iid 0.375 delta_2 + 0.625 delta_0 == G
# Objective is to estimate the density of G.
require(Rmosek)
n <- 20000
m <- 300
v <- rep(0,n)
v[sample(n,7500)] <- 2
y <- rnorm(n) + v
z <- GLmix(y, hist = TRUE, verb = 5)
plot(z,xlab = expression(mu),main = "Estimated Mixing Density")
