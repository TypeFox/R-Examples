# Test problem to compare mosek and pogs implementations
# Test problem:  y ~ N(mu_i , 1), mu_i ~ iid 0.375 delta_2 + 0.625 delta_0 == G
n <- 20000
m <- 300
v <- rep(0,n)
v[sample(n,7500)] <- 2
y <- rnorm(n) + v
lgLik <- rep(0,4)
cpu <- rep(0,4)

cpu[1] <- system.time(z <- GLmix(y, hist = TRUE, verb = 5))[3]
plot(z,xlab = expression(mu),main = "Estimated Mixing Density")
lgLik[1] <- z$logLik
leg <- as.list(1:3)
maxit <- c(10000,100000,1000000)
for(i in 1:length(maxit)){
    cpu[i+1] <- system.time(z <- GLmix(y, hist = TRUE, method = "pogs", rel_tol = 1e-8, 
	       abs_tol = 1e-8, max_iter = maxit[i]))[3]
    lgLik[i+1] <- z$logLik - lgLik[1]
    lines(z, col = i + 1)
    leg[[i]] <- paste("It = ", maxit[i], " cpu = ", 
		      round(cpu[i+1]), "L = ", round(lgLik[i+1],2))
}
legend(.5,13, leg, lty = 1, col = 2:4)
