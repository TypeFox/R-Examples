gs4 <-
function (ystar, x, idx, distr.k, sigma.k, distr.p0, mu.p0, sigma.p0) 
{
    r <- length(ystar)
    nstar <- as.numeric(table(idx))
    for (j in seq(r)) {
        id <- which(!is.na(match(idx, j)))
        xj <- x[id]
        xbar <- sum(xj)/nstar[j]
        y2star <- rk(1, distr = distr.k, mu = xbar, sigma = sigma.k/sqrt(nstar[j]))
        f.ratio <- rfystar(y2star, ystar[j], xj, distr.k = distr.k, 
            sigma.k = sigma.k, distr.p0 = distr.p0, mu.p0 = mu.p0, 
            sigma.p0 = sigma.p0)
        k.ratio <- dk(ystar[j], distr = distr.k, mu = xbar, sigma = sigma.k/sqrt(nstar[j]))/dk(y2star, 
            distr = distr.k, mu = xbar, sigma = sigma.k/sqrt(nstar[j]))
        q2 <- min(1, f.ratio * k.ratio)
        ystar[j] <- ifelse(runif(1) <= q2, y2star, ystar[j])
    }
    return(ystar)
}
