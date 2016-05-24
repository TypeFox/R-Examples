
 library(mcmc)

 set.seed(42)

 d <- 3
 witch.which <- 1 - (1 / 2)^(1 / d) * (1 / 4)^(seq(0, 5) / d)
 witch.which

 ncomp <- length(witch.which)

 neighbors <- matrix(FALSE, ncomp, ncomp)
 neighbors[row(neighbors) == col(neighbors) + 1] <- TRUE
 neighbors[row(neighbors) == col(neighbors) - 1] <- TRUE
 neighbors[row(neighbors) == col(neighbors) + 2] <- TRUE
 neighbors[row(neighbors) == col(neighbors) - 2] <- TRUE

 ludfun <- function(state) {
     stopifnot(is.numeric(state))
     stopifnot(length(state) == d + 1)
     icomp <- state[1]
     stopifnot(icomp == as.integer(icomp))
     stopifnot(1 <= icomp && icomp <= ncomp)
     theta <- state[-1]
     if (any(abs(theta) > 1.0)) return(-Inf)
     bnd <- witch.which[icomp]
     if(bnd >= 1.0)
         stop(sprintf("witch.which[%d] >= 1.0", icomp))
     if(bnd <= 0.0)
         stop(sprintf("witch.which[%d] <= 0.0", icomp))
     if (all(abs(theta) > bnd))
         return(- (d + 1) * log(2) - d * log(1 - bnd))
     return(- (d + 1) * log(2) - log1p(- (1 - bnd)^d))
 }

 initial <- c(1, rep(0, d))

 out <- temper(ludfun, initial = initial, neighbors = neighbors,
     nbatch = 50, blen = 13, nspac = 7, scale = 0.3456789)

 names(out)

 out$acceptx

 out$accepti

 colMeans(out$ibatch)

 ### check that have prob 1 / 2 for corners

 outfun <- function(state) {
     stopifnot(is.numeric(state))
     icomp <- state[1]
     stopifnot(icomp == as.integer(icomp))
     stopifnot(1 <= icomp && icomp <= length(witch.which))
     theta <- state[-1]
     foo <- all(abs(theta) > witch.which[icomp])
     bar <- rep(0, length(witch.which))
     baz <- rep(0, length(witch.which))
     bar[icomp] <- as.numeric(foo)
     baz[icomp] <- 1
     return(c(bar, baz))
 }

 out <- temper(out, blen = 103, outfun = outfun, debug = TRUE)

 eta.batch <- out$batch[ , seq(1, ncomp)]
 noo.batch <- out$batch[ , seq(ncomp + 1, ncomp + ncomp)]
 eta <- colMeans(eta.batch)
 noo <- colMeans(noo.batch)
 mu <- eta / noo
 eta
 noo
 mu

 eta.batch.rel <- sweep(eta.batch, 2, eta, "/")
 noo.batch.rel <- sweep(noo.batch, 2, noo, "/")
 mu.batch.rel <- eta.batch.rel - noo.batch.rel

 mu.mcse.rel <- apply(mu.batch.rel, 2, sd) / sqrt(out$nbatch)
 mu.mcse.rel

 foo <- cbind(mu, mu * mu.mcse.rel)
 colnames(foo) <- c("means", "MCSE")
 foo

 ### check decision about within-component or jump/swap

 identical(out$unif.which < 0.5, out$which)

 identical(out$which, out$proposal[ , 1] == out$state[ , 1])

 ### check hastings ratio calculated correctly

 n <- apply(neighbors, 1, sum)
 i <- out$state[ , 1]
 istar <- out$proposal[ , 1]
 foo <- apply(out$state, 1, ludfun)
 bar <- apply(out$proposal, 1, ludfun)
 my.log.hastings <- bar - foo - log(n[istar]) + log(n[i])
 all.equal(my.log.hastings, out$log.hastings)

