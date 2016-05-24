
as.mcmc.list.bugs <- function(x, ...)
{
    if(!inherits(x, "bugs")) 
        stop("Method as.mcmc.list.bugs() is only intended for bugs objects.")
    if(dim(x$sims.array)[2] != x$n.chains)
        stop("Inconsistancy in bug object regarding the number of chains.")
    mclis <- vector("list", x$n.chains)
    strt <- x$n.burnin + 1
    end <- x$n.iter
    ord <- order(dimnames(x$sims.array)[[3]])
    for(i in 1:x$n.chains) {
        tmp1 <- x$sims.array[,i,ord]
        mclis[[i]] <- mcmc(tmp1, start=strt, end=end, thin=x$n.thin)
    }
    as.mcmc.list(mclis)
} 