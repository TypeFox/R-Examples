sm.pca <- function(x, Y, h, cor = TRUE, nperm = 100, pc = 1, ...) {

   opt <- sm.options(list(...))
   replace.na(opt, test,    TRUE)
   replace.na(opt, band,    TRUE)
   replace.na(opt, display, c("eigenvalues", "eigenvectors"))
   replace.na(opt, col,     "red")
   
   # If x is a previously computed sm.pca object then plot as required
   
   if (is.list(x) & ("evals" %in% names(x))) {
   	
   	  evalmat <- t(x$evals)
   	  evec    <- x$evecs
      zgrid   <- x$eval.points
   	  ngrid   <- length(zgrid)
    	  if ("evals.perm" %in% names(x)) neweval <- x$evals.perm
   	  if ("evecs.perm" %in% names(x)) newevec <- x$evecs.perm
      ylim.missing <- any(is.na(opt$ylim))

      if ("eigenvalues" %in% opt$display) {
         e   <- evalmat[pc, ]
         rng <- range(e)
         if (opt$band & ("evals.perm" %in% names(x))) {
            r           <- neweval[pc, , ]
            neweig1     <- apply(r, 1, function(x) quantile(x, prob = 0.025))
            neweig2     <- apply(r, 1, function(x) quantile(x, prob = 0.975))
            rng         <- range(rng, neweig1, neweig2)
            mtch        <- match("band", names(x))
            x$band      <- cbind(q1 = neweig2, q2 = neweig1)
         }
         if (ylim.missing) opt$ylim <- rng
         plot(zgrid, e, type = "n", xlab = opt$xlab, ylab = "Variance", ylim = opt$ylim)
         if (opt$band & ("evals.perm" %in% names(x)))
            polygon(c(zgrid, rev(zgrid)), c(neweig1, rev(neweig2)), col = opt$col.band, border = NA)
         lines(zgrid, e, lwd = opt$lwd, col = opt$col)
      }

      if ("eigenvectors" %in% opt$display) {

         e <- evec[ , pc, ]
         p <- dim(e)[1]

         if (opt$band & ("evecs.perm" %in% names(x))) {
            r <- newevec[ , pc, , ]
            for(k in 1:p) {
            for(i in 1:nperm) {
            for(j in 1:ngrid) {
               if((sign(sum(e[k, ])) == sign(r[k, j, i])) == FALSE)
                  r[k, j, i] <- -1 * r[k, j, i]
            }}}

            xmat   <- matrix(nrow = 2 * ngrid - 1, ncol = p)
            clrmat <- matrix(nrow = 2 * ngrid - 2, ncol = p)
            for (i in 1:p) {
               q    <- apply(r[i, , ], 1, function(x) quantile(x, prob = c(0.025, 0.975)))
               ind  <- as.numeric(sign((q[1, ] - e[i, ]) * (q[2, ] - e[i, ])) < 1)
               ind  <- rep(ind, each = 2)
               ind  <- ind[-c(1, 2 * ngrid)]
               clrb <- rep(255, 3)
               clri <- col2rgb(i)
               clri <- rgb(clri[1] + ind * (clrb[1] - clri[1]) * 0.9,
                           clri[2] + ind * (clrb[2] - clri[2]) * 0.9,
                           clri[3] + ind * (clrb[3] - clri[3]) * 0.9,
                           maxColorValue = 255)
               xg   <- c(rbind(e[i, ], e[i, ] + c(0.5 * diff(e[i, ]), 0)))
               xmat[ , i]   <- xg[-2 * ngrid]
               clrmat[ , i] <- clri
            }
      
            zg   <- c(rbind(zgrid, zgrid + c(0.5 * diff(zgrid), 0)))
            zg   <- zg[-2 * ngrid]
            if (ylim.missing) opt$ylim <- range(xmat)
            plot(range(zg), range(xmat), type = "n", xlab = opt$xlab, ylab = "PC loadings", ylim = opt$ylim)
            for (i in 1:p)
               segments(zg[-length(zg)], xmat[-nrow(xmat), i], zg[-1], xmat[-1, i], col = clrmat[ , i], lty = i)
            x$xgrid.plot <- zg
            x$evecs.plot <- xmat
            x$col.plot   <- clrmat
         }
         else {
            replace.na(opt, ylim, range(e))
            plot(range(zgrid), range(e), type = "n", xlab = opt$xlab, ylab = "PC loadings", ylim = opt$ylim)
            for (i in 1:p) lines(zgrid, e[i, ], col = i, lty = i)
         }
 
      }
   
      if (opt$test & ("evals.perm" %in% names(x))) {
         n    <- length(x$x)
         p    <- ncol(x$Y)
         S    <- sm.weight(x$x, x$x, h = x$h, options = list(poly.index = 1))
         mht  <- S %*% x$Y
         S    <- x$Y - mht
         # S    <- sweep(x$Y, 2, apply(x$Y, 2, mean))
         sdv  <- if (cor) apply(x$Y, 2, sd) else rep(1, p)
         S    <- S / matrix(rep(sdv, each = n), ncol = p)
         S    <- t(S) %*% S / n
         pc0  <- eigen(S, symmetric = TRUE)
         
         e    <- evalmat[pc, ]
         r    <- neweval[pc, , ]         
         tst1 <- sum(abs(e - pc0$values[pc]))
         ref1 <- apply(abs(r - pc0$values[pc]), 2, sum)
         pv1  <- length(which(ref1 > tst1)) / length(ref1)
         
         e0   <- pc0$vectors[ , pc]
         e    <- evec[ , pc, ]
         r    <- newevec[ , pc, , ]
         tst2 <- sum(1 - apply(sweep(e, 1, e0, "*"), 2, sum)^2)
         ref2 <- apply(1 - apply(sweep(r, 1, e0, "*"), 2:3, sum)^2, 2, sum)
         pv2  <- length(which(ref2 > tst2)) / length(ref2)

         x$p.values  <- pv1
         x$p.vectors <- pv2
         if (opt$verbose > 0) {
            cat("Eigenvalues:  p =", round(pv1, 3), "\n")
            cat("Eigenvectors: p =", round(pv2, 3), "\n")
         }
      }
   
      return(invisible(x))
   }
   
   #  Compute the information needed for plotting
 
   if (!is.vector(x)) stop("x should be a vector.")
   if (!is.matrix(Y)) stop("Y should be a matrix.")
   if (!(length(x) == nrow(Y))) stop("the length of x should match the number of rows of Y.")
   
   x.name  <- deparse(substitute(x))
   y.name  <- deparse(substitute(Y))
   z       <- x

   ind     <- apply(cbind(x, Y), 1, function(x) any(is.na(x)))
   z       <- z[!ind]
   Y       <- Y[!ind, ]
   
   if(missing(h)) h <- h.select(x, Y[ , 1], ...)

   replace.na(opt, ngrid, 25)
   replace.na(opt, eval.points, seq(min(z), max(z), length = opt$ngrid))
   replace.na(opt, xlab,  x.name)
   zgrid <- opt$eval.points
   ngrid <- opt$ngrid

   n          <- nrow(Y)
   p          <- ncol(Y)
   evec       <- array(0, c(p, p, ngrid))
   evalmat    <- matrix(0, nrow = p, ncol = ngrid)
   varmat     <- matrix(0, nrow = p, ncol = ngrid)
   mhat       <- matrix(0, nrow = p, ncol = ngrid)
   indmat     <- matrix(0, nrow = p, ncol = ngrid)
   indmat[,1] <- 1:p
   max.ind    <- function(x, p)  which(abs(x) == max(abs(x)))

   # Vtot <- matrix(0, p, p)

   for (i in 1:ngrid) {

      S         <- sm.weight(z, zgrid[i], h = h, options = list(poly.index = 1))
      mhat[, i] <- as.vector(S %*% Y)

      D         <- Y - matrix(rep(mhat[, i], n), ncol = p, byrow = TRUE)
      sdvec     <- if (cor) apply(Y, 2, sd) else rep(1, p)
      D         <- D / matrix(rep(sdvec, each = n), ncol = p)
      S         <- sm.weight(z, zgrid[i], h = h, options = list(poly.index = 1))
      V         <- t(D) %*% diag(as.vector(S)) %*% D
	
      # Vtot      <- Vtot + length(z) * V

      pca       <- eigen(V, symmetric = TRUE)

      if (i > 1) {
        d          <- crossprod(pca$vectors, old.vectors)
        indmat[,i] <- apply(d, 1, max.ind, p = p)
        indmat[,i] <- indmat[ , i - 1][indmat[ , i]]
        dmax       <- d[cbind(1:p, indmat[ , i])]
        for (j in (1:p)[dmax < 0]) 
           pca$vectors[ , j] <- -pca$vectors[ , j]
      }
      
      pca$values[pca$values < 0] <- 0
      pcvar        <- pca$values / sum(pca$values)
      ord          <- order(indmat[ , i])
      varmat[, i]  <- pcvar[ord]
      evalmat[, i] <- pca$values[ord]
      evec[ , , i] <- pca$vectors[, ord]

      old.vectors  <- pca$vectors
      ind.old      <- indmat[,i]
   }
     
   result <- list(xgrid = zgrid, evecs = evec, evals = t(evalmat), mhat = t(mhat),
                  var.explained = t(varmat), eval.points = opt$eval.points, xlab = opt$xlab,
                  h = h, x = x, Y = Y, nperm = nperm, cor = cor)

   if (opt$test | opt$band) {
      newevec <- array(NA, c(ncol(Y), ncol(Y), ngrid, nperm))
      neweval <- array(NA, c(ncol(Y), ngrid, nperm))
      # Vtot    <- array(NA, c(dim(Y)[2], dim(Y)[2], nperm))
      znew    <- replicate(nperm, sample(z, length(z)), simplify = TRUE)

      spc <- function(x){
         spc2    <- sm.pca(x, Y, h, cor = cor, eval.points = opt$eval.points,
                                 test = FALSE, band = FALSE, display = "none")	
         # Vtot    <- spc2$Vtot
         newevec <- spc2$evecs
         neweval <- t(spc2$evals)
         list(newevec = newevec, neweval = neweval)
      }

      res <- apply(znew, 2, spc)

      for(i in 1:nperm){
         # Vtot[ , , i]    <- res[[i]]$Vtot
         newevec[ , , , i] <- res[[i]]$newevec
         neweval[ , , i]   <- res[[i]]$neweval
      }
      
   result$evecs.perm <- newevec
   result$evals.perm <- neweval
   }
   
   if (!("none" %in% opt$display)) result <- sm.pca(result, display = opt$display, pc = pc)


   invisible(result)
}

