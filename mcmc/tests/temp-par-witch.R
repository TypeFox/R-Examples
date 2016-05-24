
 library(mcmc)

 options(digits=4) # avoid rounding differences

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

 thetas <- matrix(0, ncomp, d)
 out <- temper(ludfun, initial = thetas, neighbors = neighbors, nbatch = 50,
     blen = 13, nspac = 7, scale = 0.3456789, parallel = TRUE)

 names(out)

 out$acceptx

 out$accepti

 ### check that have prob 1 / 2 for corners

 outfun <- function(state) {
     stopifnot(is.matrix(state))
     ncomp <- nrow(state)
     d <- ncol(state)
     foo <- sweep(abs(state), 1, witch.which)
     bar <- apply(foo > 0, 1, all) 
     return(as.numeric(bar))
 }

 out <- temper(out, outfun = outfun)

 colMeans(out$batch)
 apply(out$batch, 2, sd) / sqrt(out$nbatch)

 ### try again

 out <- temper(out, blen = 103, outfun = outfun)

 foo <- cbind(colMeans(out$batch),
     apply(out$batch, 2, sd) / sqrt(out$nbatch))
 colnames(foo) <- c("means", "MCSE")
 foo

