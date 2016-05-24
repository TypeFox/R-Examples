
 library(potts)

 set.seed(42)

 ncolor <- as.integer(4)
 alpha <- rnorm(ncolor) * 0.01
 beta <- log(1 + sqrt(ncolor))
 theta <- c(alpha, beta)

 nrow <- 25
 ncol <- 20
 x <- matrix(1, nrow = nrow, ncol = ncol)
 foo <- packPotts(x, ncolor)

 outfun <- function(tt) {
     qux <- outer(tt, tt)
     c(tt, qux[lower.tri(qux, diag = TRUE)])
 }

 outfun(c(485, 2, 9, 4, 954))

 out <- potts(foo, theta, nbatch = 5, blen = 3, nspac = 2, debug = TRUE,
     outfun = outfun)
 names(out)

 niter <- out$nbatch * out$blen * out$nspac

 .Random.seed <- out$initial.seed
 out.too <- potts(foo, theta, nbatch = niter)

 tt <- out.too$batch
 ttaug <- t(apply(tt, 1, outfun))
 identical(tt, ttaug[ , 1:ncol(tt)])
 nout <- ncol(ttaug)
 ncol(out$batch) == nout

 ##### check batch means #####

 foo <- ttaug[seq(1, niter) %% out$nspac == 0, ]
 foo <- array(as.vector(foo), c(out$blen, out$nbatch, nout))
 foo <- apply(foo, c(2, 3), mean)
 identical(foo, out$batch)

