Haz3 <- function (d, S, gam, lam, fvec) 
{
    gam12 <- as.vector(gam[, 1] - gam[, 2])
    temp <- gam12 * S
    gweight <- (1 - d * temp)/(gam[, 2] + temp)
    n <- length(d)
    Hw <- rep(0, n)
    ttemp <- cumsumsurv(gweight) + lam*fvec    
	   ###  rev(cumsum(rev(gweight))) + lam*fvec  3/2015 changed by MZ
    Hw <- d/ttemp
###################### Change Nov. 2013.
#    for (k in 1:n) {
#        if (d[k] == 1) 
#            Hw[k] <- 1/(sum(gweight[k:n]) + lam * fvec[k])
#    }
###################### get rid of for loop
    list(Hazw = Hw, Su = exp(-cumsum(Hw)))
}
