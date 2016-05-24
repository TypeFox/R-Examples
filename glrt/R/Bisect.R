Bisect <-
function (tA, pvec, ndir, Meps, tolbis = 1e-07) 
{
    etainv <- 1/(tA %*% pvec)
    bot <- 0
    top <- 1
    mult <- tA %*% ndir
    dbot <- sum(etainv * mult)
    ptop <- rescaleP(pvec + top * ndir, Meps)
    pbot <- pvec
    dtop <- sum(mult/(tA %*% ptop))
    done <- FALSE
    while (!done) {
        if (sign(dbot) * sign(dtop) > 0 || top - bot < tolbis) {
            ltop <- sum(log(tA %*% ptop))
            lbot <- sum(log(tA %*% pbot))
            if (lbot > ltop) 
                pnew <- rescaleP(pvec + bot * ndir, Meps)
            else pnew <- rescaleP(pvec + top * ndir, Meps)
            done <- TRUE
        }
        else {
            mid <- (bot + top)/2
            pmid <- rescaleP(pvec + mid * ndir, Meps)
            dmid <- sum(mult/(tA %*% pmid))
            if (dmid * dtop < 0) {
                bot <- mid
                dbot <- dmid
                pbot <- pmid
            }
            else {
                top <- mid
                dtop <- dmid
                ptop <- pmid
            }
        }
    }
    pnew
}