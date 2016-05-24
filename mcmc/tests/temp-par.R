
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
 betas <- NULL
 for (i in 1:nrow(models)) {
     inies <- as.logical(models[i, ])
     foo <- exes[inies]
     bar <- paste("y ~", paste(foo, collapse = " + "))
     if (! any(inies)) bar <- paste(bar, "1")
     baz <- glm(as.formula(bar), family = binomial)
     beta <- rep(0, 4)
     beta[c(TRUE, inies)] <- baz$coef
     betas <- rbind(betas, beta)
 }

 neighbors <- matrix(FALSE, nrow(models), nrow(models))
 for (i in 1:nrow(neighbors)) {
     for (j in 1:ncol(neighbors)) {
         foo <- models[i, ]
         bar <- models[j, ]
         if (sum(foo != bar) == 1) neighbors[i, j] <- TRUE
     }
 }

 ludfun <- function(state, ...) {
     stopifnot(is.numeric(state))
     stopifnot(length(state) == ncol(models) + 2)
     stopifnot(length(state) == ncol(models) + 2)
     icomp <- state[1]
     stopifnot(icomp == as.integer(icomp))
     stopifnot(1 <= icomp && icomp <= nrow(models))
     beta <- state[-1]
     inies <- c(TRUE, as.logical(models[icomp, ]))
     beta.logl <- beta
     beta.logl[! inies] <- 0
     eta <- as.numeric(modmat %*% beta.logl)
     logp <- ifelse(eta < 0, eta - log1p(exp(eta)), - log1p(exp(- eta)))
     logq <- ifelse(eta < 0, - log1p(exp(eta)), - eta - log1p(exp(- eta)))
     logl <- sum(logp[y == 1]) + sum(logq[y == 0])
     val <- logl - sum(beta^2) / 2
     return(val)
 }

 ludval <- NULL
 for (i in 1:nrow(models)) ludval <- c(ludval, ludfun(c(i, betas[i, ])))
 all(is.finite(ludval))


 out <- temper(ludfun, initial = betas, neighbors = neighbors, nbatch = 20,
     blen = 10, nspac = 5, scale = 0.56789, parallel = TRUE, debug = TRUE)

 names(out)

 ### check decision about within-component or jump/swap

 identical(out$unif.which < 0.5, out$which)

 identical(out$which, out$proposal[ , 1] == out$coproposal[ , 1])

 ### check proposal and coproposal are actually current state or part thereof

 prop <- out$proposal
 coprop <- out$coproposal
 prop.i <- prop[ , 1]
 coprop.i <- coprop[ , 1]
 alt.prop <- prop
 alt.coprop <- coprop
 for (i in 1:nrow(prop)) {
     alt.prop[i, ] <- c(prop.i[i], out$state[i, prop.i[i], ])
     alt.coprop[i, ] <- c(coprop.i[i], out$state[i, coprop.i[i], ])
 }
 identical(coprop, alt.coprop)
 identical(prop[! out$which, ], alt.prop[! out$which, ])
 identical(prop[out$which, 1], alt.prop[out$which, 1])

 ### check hastings ratio calculated correctly

 foo <- apply(prop, 1, ludfun)
 fooco <- apply(coprop, 1, ludfun)
 prop[ , 1] <- out$coproposal[ , 1]
 coprop[ , 1] <- out$proposal[ , 1]
 foo.swap <- apply(prop, 1, ludfun)
 fooco.swap <- apply(coprop, 1, ludfun)
 log.haste <- ifelse(out$which, foo - fooco,
     foo.swap + fooco.swap - foo - fooco)
 all.equal(log.haste, out$log.hastings)

 ### check hastings rejection decided correctly

 identical(out$log.hastings >= 0, is.na(out$unif.hastings))
 all(out$log.hastings < 0 | out$acceptd)
 identical(out$acceptd,
     out$log.hastings >= 0 | out$unif.hastings < exp(out$log.hastings))

 ### check acceptance carried out or not (according to decision) correctly

 before <- out$state
 after <- before
 after[- dim(after)[1], , ] <- before[-1, , ]
 after[dim(after)[1], , ] <- out$final
 my.after <- before
 for (i in 1:length(out$acceptd)) {
     if (out$acceptd[i]) {
         if (out$which[i]) {
             j <- out$proposal[i, 1]
             my.after[i, j, ] <- out$proposal[i, -1]
         } else {
             j <- out$proposal[i, 1]
             k <- out$coproposal[i, 1]
             my.after[i, j, ] <- out$coproposal[i, -1]
             my.after[i, k, ] <- out$proposal[i, -1]
         }
     }
 }
 identical(after, my.after)

 ### check within-component proposal

 my.coproposal.within <- out$coproposal[out$which, ]
 proposal.within <- out$proposal[out$which, ]
 my.z <- out$norm[out$which, ]
 my.proposal.within <- my.coproposal.within
 my.proposal.within[ , -1] <- my.coproposal.within[ , -1] + out$scale * my.z
 all.equal(proposal.within, my.proposal.within)

 my.unif.choose <- out$unif.choose[out$which, 1]
 my.i <- floor(nrow(models) * my.unif.choose) + 1
 all(1 <= my.i & my.i <= nrow(models))
 identical(my.i, my.coproposal.within[ , 1])

 ### check swap proposal

 coproposal.swap <- out$coproposal[! out$which, ]
 proposal.swap <- out$proposal[! out$which, ]
 unif.choose.swap <- out$unif.choose[! out$which, ]
 my.i <- floor(nrow(models) * unif.choose.swap[ , 1]) + 1
 nneighbors <- apply(out$neighbors, 1, sum)
 my.nneighbors <- nneighbors[my.i]
 my.k <- floor(my.nneighbors * unif.choose.swap[ , 2]) + 1
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

 nx <- ncol(out$initial)
 niter <- out$nbatch * out$blen * out$nspac
 my.norm <- matrix(NA, nrow = nrow(out$norm), ncol = ncol(out$norm))
 my.unif.which <- rep(NA, niter)
 my.unif.hastings <- rep(NA, niter)
 my.unif.choose <- matrix(NA, niter, 2)
 for (iiter in 1:niter) {
     my.unif.which[iiter] <- runif(1)
     if (out$which[iiter]) {
         my.unif.choose[iiter, 1] <- runif(1)
         my.norm[iiter, ] <- rnorm(nx)
         if (out$log.hastings[iiter] < 0) my.unif.hastings[iiter] <- runif(1) 
     } else {
         my.unif.choose[iiter, ] <- runif(2)
         if (out$log.hastings[iiter] < 0) my.unif.hastings[iiter] <- runif(1) 
     }
 }
 identical(my.norm, out$norm)
 identical(my.unif.which, out$unif.which)
 identical(my.unif.hastings, out$unif.hastings)
 identical(my.unif.choose, out$unif.choose)

 .Random.seed <- save.Random.seed

 ### check batch means

 foo <- after[seq(1, niter) %% out$nspac == 0, , ]
 foo <- array(as.vector(foo), dim = c(out$blen, out$nbatch, dim(foo)[2:3]))
 foo <- apply(foo, c(2, 3, 4), mean)
 all.equal(foo, out$batch)

 ### check acceptance rates

 accept.within <- out$acceptd[out$which]
 my.i.within <- out$coproposal[out$which, 1]
 my.acceptx <- as.vector(sapply(split(accept.within, my.i.within), mean))
 identical(my.acceptx, out$acceptx)

 accept.swap <- out$acceptd[! out$which]
 my.i.swap <- out$coproposal[! out$which, 1]
 my.j.swap <- out$proposal[! out$which, 1]
 nmodel <- nrow(out$neighbors)
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

 out <- temper(out, scale = newscale)

 my.coproposal.within <- out$coproposal[out$which, ]
 proposal.within <- out$proposal[out$which, ]
 my.z <- out$norm[out$which, ]
 my.proposal.within <- my.coproposal.within
 my.proposal.within[ , -1] <- my.coproposal.within[ , -1] +
     sweep(my.z, 2, out$scale, "*")
 all.equal(proposal.within, my.proposal.within)

 ### check scale matrix

 matscale <- matrix(rnorm(nx * nx, 0.0, 0.1), nx, nx)
 diag(matscale) <- 0.56789

 out <- temper(out, scale = matscale)

 my.coproposal.within <- out$coproposal[out$which, ]
 proposal.within <- out$proposal[out$which, ]
 my.z <- out$norm[out$which, ]
 my.proposal.within <- my.coproposal.within
 my.proposal.within[ , -1] <- my.coproposal.within[ , -1] +
     my.z %*% t(out$scale)
 all.equal(proposal.within, my.proposal.within)

 ### check scale list

 lisztscale <- list(0.56789, newscale, matscale, matscale, newscale,
     0.98765, 0.98765, newscale)

 out <- temper(out, scale = lisztscale)

 my.coproposal.within <- out$coproposal[out$which, ]
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

 outfun <- function(state, icomp, ...) {
     stopifnot(is.matrix(state))
     stopifnot(is.numeric(state))
     nx <- ncol(betas)
     ncomp <- nrow(betas)
     stopifnot(ncol(state) == nx)
     stopifnot(nrow(state) == ncomp)
     stopifnot(1 <= icomp && icomp <= ncomp)
     foo <- state[icomp, ]
     bar <- foo^2
     return(c(foo, bar))
 }

 out <- temper(out, outfun = outfun, icomp = 4)

 before <- out$state
 after <- before
 after[- dim(after)[1], , ] <- before[-1, , ]
 after[dim(after)[1], , ] <- out$final
 outies <- apply(after, 1, outfun, icomp = 4)
 outies <- t(outies)

 foo <- outies[seq(1, niter) %% out$nspac == 0, ]
 foo <- array(as.vector(foo), dim = c(out$blen, out$nbatch, dim(foo)[2]))
 foo <- apply(foo, c(2, 3), mean)
 all.equal(foo, out$batch)

