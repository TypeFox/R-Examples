as.mcmc.bugs <- function(x, ...){
    n.chains <- x$n.chains
    sims <- x$sims.array
    n.thin <- x$n.thin
    if (n.chains==1) return(coda::mcmc(sims[, 1, ], thin=n.thin))
    out <- vector("list", length=n.chains)
    for (i in seq(n.chains)) out[[i]] <- mcmc(sims[, i, ], thin=n.thin)
    out <- mcmc.list(out)
    varnames(out) <- dimnames(sims)[[3]]
    return(out)
}

