# copy from coda

as.mcmc.rjags <- function (x, ...) {
  x <- x$BUGSoutput
  mclist <- vector("list", x$n.chains)
  mclis <- vector("list", x$n.chains)
  strt <- x$n.burnin + 1
  end <- x$n.iter
  ord <- order(dimnames(x$sims.array)[[3]])
  for (i in 1:x$n.chains) {
    tmp1 <- x$sims.array[, i, ord]
    mclis[[i]] <- mcmc(tmp1, start = strt, end = end, thin = x$n.thin)
  }
  as.mcmc.list(mclis)
}




 # if (x$n.chains > 1) {
#    z <- list()
#    for (i in 1:x$n.chains) {
#      z[[i]] <- mcmc(x$sims.array[, i, ], start = 1, thin = x$n.thin)
#    }
#  class(z) <- "mcmc.list"
#  }
#  else {
#    z <- mcmc(x$sims.matrix, start = 1, thin = x$n.thin)
#  }
#  return(z)
