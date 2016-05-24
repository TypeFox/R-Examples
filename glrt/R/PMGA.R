PMGA <- 
function (est, ww = rep(1, length(est))) 
{
    curm <- length(est)
    poolnum <- rep(1, curm)
    passes <- 0
    iso <- FALSE
    while (!iso) {
        iso <- TRUE
        poolind <- 1
        curind <- 1
        while (curind <= curm) {
            groupstart <- curind
            while (curind < curm && est[curind + 1] < est[curind]) curind <- curind + 
                1
            iso <- poolind == curind
            est[poolind] <- sum(ww[groupstart:curind] * est[groupstart:curind])/sum(ww[groupstart:curind])
            ww[poolind] <- sum(ww[groupstart:curind])
            poolnum[poolind] <- sum(poolnum[groupstart:curind])
            poolind <- poolind + 1
            curind <- curind + 1
        }
        curm <- poolind - 1
        passes <- passes + 1
    }
    return(list(est = est[1:curm], ww = ww[1:curm], poolnum = poolnum[1:curm], 
        passes = passes))
}