bvn.plot <- function (mu = c(0, 0), cv = 0, vr = c(1, 1), res = 0.3, xlab = expression(y[1]), ylab = expression(y[2]), 
        zlab = expression(paste("f(", y[1], ",", y[2], ")")),...) 
{
   # require(lattice::wireframe) || stop("function requires package lattice")
  
    mv <- max(vr)
    g <- expand.grid(x = seq(mu[1] - 3 * sqrt(mv), mu[1] + 3 * 
        sqrt(mv), res), y = seq(mu[2] - 3 * sqrt(mv), mu[2] + 
        3 * sqrt(mv), res))
    g$z <- dmvnorm(cbind(g[, 1], g[, 2]), mu, sigma = matrix(nrow = 2, 
        ncol = 2, data = c(vr[1], cv, cv, vr[2])))
    wireframe(z ~ x * y, data = g, xlab = xlab, ylab = ylab, zlab = zlab, 
        ...)
}