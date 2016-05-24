
 library(mcmc)

 set.seed(42)

 data(foo)
 attach(foo)

 out <- glm(y ~ x1 + x2 + x3, family = binomial, x = TRUE)
 summary(out)

 modmat <- out$x

 models <- cbind(rep(0:1, each = 4), rep(rep(0:1, times = 2), each = 2),
               rep(0:1, times = 4))

 exes <- paste("x", 1:3, sep = "")
 models[nrow(models), ]
 beta.initial <- c(nrow(models), out$coefficients)

 neighbors <- matrix(FALSE, nrow(models), nrow(models))
 for (i in 1:nrow(neighbors)) {
     for (j in 1:ncol(neighbors)) {
         foo <- models[i, ]
         bar <- models[j, ]
         if (sum(foo != bar) == 1) neighbors[i, j] <- TRUE
     }
 }
 neighbors

 ludfun <- function(state, log.pseudo.prior, ...) {
     stopifnot(is.numeric(state))
     stopifnot(length(state) == ncol(models) + 2)
     icomp <- state[1]
     stopifnot(icomp == as.integer(icomp))
     stopifnot(1 <= icomp && icomp <= nrow(models))
     stopifnot(is.numeric(log.pseudo.prior))
     stopifnot(length(log.pseudo.prior) == nrow(models))
     beta <- state[-1]
     inies <- c(TRUE, as.logical(models[icomp, ]))
     beta.logl <- beta
     beta.logl[! inies] <- 0
     eta <- as.numeric(modmat %*% beta.logl)
     logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
     logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
     logl <- sum(logp[y == 1]) + sum(logq[y == 0])
     val <- logl - sum(beta^2) / 2 + log.pseudo.prior[icomp]
     return(val)
 }

 qux <- c(25.01, 5.875, 9.028, 0.6959, 11.73,  2.367, 5.864, 0.0)

 out <- temper(ludfun, initial = beta.initial, neighbors = neighbors,
     nbatch = 25, blen = 20, nspac = 5, scale = 0.56789, debug = TRUE,
     log.pseudo.prior = qux)

 names(out)

 apply(out$ibatch, 2, mean)

 ### check decision about within-component or jump/swap

 identical(out$unif.which < 0.5, out$which)

 identical(out$which, out$proposal[ , 1] == out$state[ , 1])

 ### check hastings ratio calculated correctly

 foo <- apply(out$state, 1, ludfun, log.pseudo.prior = qux)
 bar <- apply(out$proposal, 1, ludfun, log.pseudo.prior = qux)
 all.equal(bar - foo, out$log.hastings)

 ### check hastings rejection decided correctly

 identical(out$log.hastings >= 0, is.na(out$unif.hastings))
 all(out$log.hastings < 0 | out$acceptd)
 identical(out$acceptd,
     out$log.hastings >= 0 | out$unif.hastings < exp(out$log.hastings))

 ### check acceptance carried out or not (according to decision) correctly

 before <- out$state
 after <- before
 after[- dim(after)[1], ] <- before[-1, ]
 after[dim(after)[1], ] <- out$final
 my.after <- before
 my.after[out$acceptd, ] <- out$proposal[out$acceptd, ]
 identical(after, my.after)

 ### check within-component proposal

 my.coproposal.within <- out$state[out$which, ]
 proposal.within <- out$proposal[out$which, ]
 my.z <- out$norm[out$which, ]
 my.proposal.within <- my.coproposal.within
 my.proposal.within[ , -1] <- my.coproposal.within[ , -1] + out$scale * my.z
 all.equal(proposal.within, my.proposal.within)

 ### check swap proposal

 coproposal.swap <- out$state[! out$which, ]
 proposal.swap <- out$proposal[! out$which, ]
 unif.choose.swap <- out$unif.choose[! out$which]
 my.i <- coproposal.swap[ , 1]
 nneighbors <- apply(out$neighbors, 1, sum)
 my.nneighbors <- nneighbors[my.i]
 my.k <- floor(my.nneighbors * unif.choose.swap) + 1
 my.j <- my.k
 foo <- seq(1, ncol(out$neighbors))
 for (i in seq(along = my.j)) {
     my.j[i] <- (foo[out$neighbors[my.i[i], ]])[my.k[i]]
 }
 identical(coproposal.swap[ , 1], my.i)
 identical(proposal.swap[ , 1], my.j)

 ### check standard normal and uniform random numbers are as purported

 save.Random.seed <- .Random.seed
 .Random.seed <- out$initial.seed

 nx <- length(out$initial) - 1
 niter <- out$nbatch * out$blen * out$nspac
 my.norm <- matrix(NA, nrow = nrow(out$norm), ncol = ncol(out$norm))
 my.unif.which <- rep(NA, niter)
 my.unif.hastings <- rep(NA, niter)
 my.unif.choose <- rep(NA, niter)
 for (iiter in 1:niter) {
     my.unif.which[iiter] <- runif(1)
     if (out$which[iiter]) {
         my.norm[iiter, ] <- rnorm(nx)
         if (out$log.hastings[iiter] < 0) my.unif.hastings[iiter] <- runif(1) 
     } else {
         my.unif.choose[iiter] <- runif(1)
         if (out$log.hastings[iiter] < 0) my.unif.hastings[iiter] <- runif(1) 
     }
 }
 identical(my.norm, out$norm)
 identical(my.unif.which, out$unif.which)
 identical(my.unif.hastings, out$unif.hastings)
 identical(my.unif.choose, out$unif.choose)

 .Random.seed <- save.Random.seed

 ### check batch means

 my.xstate <- after[ , -1]
 foo <- my.xstate[seq(1, niter) %% out$nspac == 0, ]
 foo <- array(as.vector(foo), dim = c(out$blen, out$nbatch, dim(foo)[2]))
 foo <- apply(foo, c(2, 3), mean)
 all.equal(foo, out$batch)

 ### check ibatch means

 my.istate <- after[ , 1]
 my.istate.matrix <- matrix(0, length(my.istate), nrow(models))
 for (i in 1:nrow(my.istate.matrix))
     my.istate.matrix[i, my.istate[i]] <- 1
 foo <- my.istate.matrix[seq(1, niter) %% out$nspac == 0, ]
 foo <- array(as.vector(foo), dim = c(out$blen, out$nbatch, dim(foo)[2]))
 foo <- apply(foo, c(2, 3), mean)
 all.equal(foo, out$ibatch)

 ### check acceptance rates

 nmodel <- nrow(out$neighbors)

 accept.within <- out$acceptd[out$which]
 my.i.within <- out$state[out$which, 1]
 my.i.within.accept <- my.i.within[accept.within]
 my.acceptx.numer <- tabulate(my.i.within.accept, nbins = nmodel)
 my.acceptx.denom <- tabulate(my.i.within, nbins = nmodel)
 my.acceptx <- my.acceptx.numer / my.acceptx.denom
 identical(my.acceptx, out$acceptx)

 accept.swap <- out$acceptd[! out$which]
 my.i.swap <- out$state[! out$which, 1]
 my.j.swap <- out$proposal[! out$which, 1]
 my.accepti <- matrix(NA, nmodel, nmodel)
 for (i in 1:nmodel) {
     for (j in 1:nmodel) {
         if (out$neighbors[i, j]) {
             my.accepti[i, j] <-
                 mean(accept.swap[my.i.swap == i & my.j.swap == j])
         }
     }
 }
 identical(my.accepti, out$accepti)

 ### check scale vector

 nx <- ncol(models) + 1
 newscale <- rnorm(nx, 0.5, 0.1)

 out <- temper(out, scale = newscale, log.pseudo.prior = qux)

 my.coproposal.within <- out$state[out$which, ]
 proposal.within <- out$proposal[out$which, ]
 my.z <- out$norm[out$which, ]
 my.proposal.within <- my.coproposal.within
 my.proposal.within[ , -1] <- my.coproposal.within[ , -1] +
     sweep(my.z, 2, out$scale, "*")
 all.equal(proposal.within, my.proposal.within)

 ### check scale matrix

 matscale <- matrix(rnorm(nx * nx, 0.0, 0.1), nx, nx)
 diag(matscale) <- 0.56789

 out <- temper(out, scale = matscale, log.pseudo.prior = qux)

 my.coproposal.within <- out$state[out$which, ]
 proposal.within <- out$proposal[out$which, ]
 my.z <- out$norm[out$which, ]
 my.proposal.within <- my.coproposal.within
 my.proposal.within[ , -1] <- my.coproposal.within[ , -1] +
     my.z %*% t(out$scale)
 all.equal(proposal.within, my.proposal.within)

 ### check scale list

 lisztscale <- list(0.56789, newscale, matscale, matscale, newscale,
     0.98765, 0.98765, newscale)

 out <- temper(out, scale = lisztscale, log.pseudo.prior = qux)

 my.coproposal.within <- out$state[out$which, ]
 proposal.within <- out$proposal[out$which, ]
 my.z <- out$norm[out$which, ]
 my.proposal.within <- my.coproposal.within
 for (iiter in 1:nrow(my.z)) {
     my.i <- my.coproposal.within[iiter, 1]
     my.scale <- out$scale[[my.i]]
     if (is.matrix(my.scale)) {
         my.proposal.within[iiter, -1] <- my.coproposal.within[iiter, -1] +
             my.z[iiter, , drop = FALSE] %*% t(my.scale)
     } else {
         my.proposal.within[iiter, -1] <- my.coproposal.within[iiter, -1] +
             my.z[iiter, ] * my.scale
     }
 }
 all.equal(proposal.within, my.proposal.within)

 ### check outfun

 outfun <- function(state, icomp) {
     stopifnot(is.matrix(state))
     stopifnot(is.numeric(state))
     nx <- ncol(initial)
     ncomp <- nrow(initial)
     stopifnot(ncol(state) == nx)
     stopifnot(nrow(state) == ncomp)
     stopifnot(1 <= icomp & icomp <= ncomp)
     foo <- state[icomp, ]
     bar <- foo^2
     return(c(foo, bar))
 }

 ncomp <- nrow(models)
 nx <- length(beta.initial) - 1

 outfun <- function(state, icomp, ...) {
     stopifnot(is.numeric(state))
     stopifnot(length(state) == nx + 1)
     istate <- state[1]
     stopifnot(istate == as.integer(istate))
     stopifnot(1 <= istate && istate <= ncomp)
     stopifnot(1 <= icomp && icomp <= ncomp)
     if (istate == icomp) {
         foo <- state[-1]
     } else {
         foo <- rep(0, nx)
     }
     bar <- foo^2
     return(c(foo, bar))
 }

 out <- temper(ludfun, initial = out$final, neighbors = neighbors,
     nbatch = 25, blen = 20, nspac = 5, scale = 0.56789, debug = TRUE,
     outfun = outfun, log.pseudo.prior = qux, icomp = 4)

 before <- out$state
 after <- before
 after[- dim(after)[1], ] <- before[-1, ]
 after[dim(after)[1], ] <- out$final
 outies <- apply(after, 1, outfun, icomp = 4)
 outies <- t(outies)

 foo <- outies[seq(1, niter) %% out$nspac == 0, ]
 foo <- array(as.vector(foo), dim = c(out$blen, out$nbatch, dim(foo)[2]))
 foo <- apply(foo, c(2, 3), mean)
 all.equal(foo, out$batch)

