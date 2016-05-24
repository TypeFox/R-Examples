gsYZstar <-
function (ystar, zstar, nstar, rstar, idx, x, delta, kappa, distr.k, 
    distr.py0, mu.py0, sigma.py0, distr.pz0, mu.pz0, sigma.pz0) 
{
    for (j in seq(rstar)) {
        flag <- 1
        while (flag == 1) {
            id <- which(!is.na(match(idx, j)))
            xj <- x[id]
            xbar <- sum(xj)/nstar[j]
            z2star <- rk(1, distr = distr.pz0, mu = zstar[j], 
                sigma = zstar[j]/sqrt(delta))
            y2star <- rk(1, distr = distr.py0, mu = xbar, sigma = kappa * 
                z2star/sqrt(nstar[j]))
            f.ratio <- rfyzstar(y2star, ystar[j], z2star, zstar[j], 
                xj, distr.k = distr.k, distr.py0 = distr.py0, 
                mu.py0 = mu.py0, sigma.py0 = sigma.py0, distr.pz0 = distr.pz0, 
                mu.pz0 = mu.pz0, sigma.pz0 = sigma.pz0)
            k.ratioNum <- dk(zstar[j], distr = distr.pz0, mu = z2star, 
                sigma = z2star/sqrt(delta))
            k.ratioDen <- dk(z2star, distr = distr.pz0, mu = zstar[j], 
                sigma = zstar[j]/sqrt(delta))
            k.ratio <- k.ratioNum/k.ratioDen
            k.ratioNum <- dk(ystar[j], distr = distr.py0, mu = xbar, 
                sigma = kappa * zstar[j]/sqrt(nstar[j]))
            k.ratioDen <- dk(y2star, distr = distr.py0, mu = xbar, 
                sigma = kappa * z2star/sqrt(nstar[j]))
            k.ratio <- k.ratio * k.ratioNum/k.ratioDen
            q2 <- min(1, f.ratio * k.ratio)
            if (is.na(q2)) {
                flag <- 1
            }
            else {
                flag <- 0
                if (runif(1) <= q2) {
                  ystar[j] <- y2star
                  zstar[j] <- z2star
                }
            }
        }
    }
    return(list(ystar = ystar, zstar = zstar))
}
