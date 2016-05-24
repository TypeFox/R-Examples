`yuleSim` <-
function(ntaxa, nsets, lambda = 0.01)
{
    rmat <- matrix(0, nsets, (ntaxa-1))
    for (i in 1:nsets)
    {
      rmat[i, ] <- rexp(ntaxa-1, rate = lambda*(2:ntaxa))

    }
    #slight problem here, because the final branching time - or the
    #time to present after the birth of the Nth lineage - has a different
    #distribution (exp(-nrt)), rather than exp.

    for (i in 1:nsets)
    {
      rmat[i,] <- rev(cumsum(rmat[i, (ntaxa-1):1]))

    }

    return(rmat)
}

