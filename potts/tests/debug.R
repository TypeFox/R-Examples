
 library(potts)
 library(pooh)

 set.seed(42)

 ncolor <- as.integer(4)
 alpha <- rnorm(ncolor) * 0.01
 beta <- log(1 + sqrt(ncolor))
 theta <- c(alpha, beta)

 nrow <- 25
 ncol <- 20
 x <- matrix(1, nrow = nrow, ncol = ncol)
 foo <- packPotts(x, ncolor)

 out <- potts(foo, theta, nbatch = 5, blen = 3, nspac = 2, debug = TRUE)
 names(out)

 identical(out$initial, foo)

 before <- out$pstate
 dim(before)

 niter <- dim(before)[1]
 niter == out$nbatch * out$blen * out$nspac

 after <- before
 after[- niter, , ] <- before[- 1, ,]
 after[niter, , ] <- unpackPotts(out$final)

 sort(unique(as.vector(before)))
 sort(unique(as.vector(after)))
 all.equal(x, before[1, , ])

 ##### calculate canonical statistics #####

 ttt <- matrix(NA, niter, ncolor + 1)
 for (icolor in 1:ncolor)
     ttt[ , icolor] <- apply(after == icolor, 1, sum)

 tstar <- rep(0, niter)
 tstar <- tstar + apply(after[ , 1, ] == after[ , nrow, ], 1, sum)
 for (i in 2:nrow)
     tstar <- tstar + apply(after[ , i, ] == after[ , i - 1, ], 1, sum)
 tstar <- tstar + apply(after[ , , 1] == after[ , , ncol], 1, sum)
 for (i in 2:ncol)
     tstar <- tstar + apply(after[ , , i] == after[ , , i - 1], 1, sum)
 ttt[ , ncolor + 1] <- tstar

 ##### check batch means #####

 foo <- ttt[seq(1, niter) %% out$nspac == 0, ]
 foo <- array(as.vector(foo), c(out$blen, out$nbatch, ncolor + 1))
 foo <- apply(foo, c(2, 3), mean)
 identical(foo, out$batch)

 ##### check bonds #####

 bprob <- (- expm1(- beta))
 all.equal(bprob, 1 - exp(- beta))

 my.hstate.possible <- array(FALSE, c(niter, nrow, ncol))
 for (i in 1:nrow) {
     ip1 <- ifelse(i < nrow, i + 1, 1)
     my.hstate.possible[ , i, ] <- before[ , i, ] == before[ , ip1, ]
 }
 storage.mode(my.hstate.possible) <- "integer"
 identical(my.hstate.possible == 1, out$hunif != -1)
 my.hstate <- out$hunif < bprob
 storage.mode(my.hstate) <- "integer"
 my.hstate <- my.hstate * my.hstate.possible
 identical(my.hstate, out$hstate)

 my.vstate.possible <- array(FALSE, c(niter, nrow, ncol))
 for (i in 1:ncol) {
     ip1 <- ifelse(i < ncol, i + 1, 1)
     my.vstate.possible[ , , i] <- before[ , , i] == before[ , , ip1]
 }
 storage.mode(my.vstate.possible) <- "integer"
 identical(my.vstate.possible == 1, out$vunif != -1)
 my.vstate <- out$vunif < bprob
 storage.mode(my.vstate) <- "integer"
 my.vstate <- my.vstate * my.vstate.possible
 identical(my.vstate, out$vstate)

 ##### check patches #####
 
 my.row <- row(my.hstate[1, , ])
 my.col <- col(my.hstate[1, , ])
 my.other.row <- my.row + 1
 my.other.row[my.other.row > nrow] <- 1
 my.other.col <- my.col + 1
 my.other.col[my.other.col > ncol] <- 1

 vertices <- paste(my.row, my.col, sep = ":")

 patch.equals <- NULL

 for (iiter in 1:niter) {

     isbond <- my.hstate[iiter, , ] == 1
     my.row.bond <- as.vector(my.row[isbond])
     my.col.bond <- as.vector(my.col[isbond])
     my.other.row.bond <- as.vector(my.other.row[isbond])
     my.from.h <- paste(my.row.bond, my.col.bond, sep = ":")
     my.to.h <- paste(my.other.row.bond, my.col.bond, sep = ":")

     isbond <- my.vstate[iiter, , ] == 1
     my.row.bond <- as.vector(my.row[isbond])
     my.col.bond <- as.vector(my.col[isbond])
     my.other.col.bond <- as.vector(my.other.col[isbond])
     my.from.v <- paste(my.row.bond, my.col.bond, sep = ":")
     my.to.v <- paste(my.row.bond, my.other.col.bond, sep = ":")

     wout <- weak(from = c(my.from.h, my.from.v), to = c(my.to.h, my.to.v),
         domain = vertices, markers = TRUE)
     blab <- as.vector(out$patch[iiter, , ])
     widx <- sort(unique(wout))
     bidx <- blab[match(widx, wout)]

     patch.equals <- c(patch.equals, identical(bidx[wout], blab))
 }

 all(patch.equals)

 ##### check colors #####

 unif.equals <- NULL
 my.after <- after

 for (iiter in 1:niter) {
     blab <- as.vector(out$patch[iiter, , ])
     blab.count <- tabulate(blab, nbins = nrow * ncol)
     punif <- out$punif[iiter, ]
     unif.equals <- c(unif.equals, identical(punif != -1, blab.count != 0))
     alpha.prod <- outer(blab.count, alpha)
     p.prod <- exp(alpha.prod)
     p.sum <- apply(p.prod, 1, sum)
     p <- sweep(p.prod, 1, p.sum, "/")
     p.cum <- apply(p, 1, cumsum)
     p.cum <- t(p.cum)
     p.foo <- sweep(p.cum, 1, out$punif[iiter, ])
     newcolor <- apply(p.foo < 0, 1, sum) + 1
     newcolor[blab.count == 0] <- NA
     blab.new <- newcolor[blab]
     my.after[iiter, , ] <- matrix(blab.new, nrow, ncol)
 }

 all(unif.equals)
 all.equal(after, my.after)

 ##### check uniform random numbers #####

 u <- c(as.vector(out$hunif), as.vector(out$vunif), as.vector(out$punif))
 u <- u[u != -1]
 ks.test(x = u, y = "punif")

 ##### check alpha = 0 #####
 
 alpha <- rep(0, ncolor)
 theta <- c(alpha, beta)

 out.too <- potts(out, param = theta, debug = FALSE)

 alpha <- rep(0.1, ncolor)
 theta <- c(alpha, beta)

 out.too.too <- potts(out, param = theta, debug = FALSE)

 all.equal(out.too$batch, out.too.too$batch)

