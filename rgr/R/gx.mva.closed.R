gx.mva.closed <-
function(xx, main = deparse(substitute(xx)))
{
     # Special version of mva for use with closed, compositional, data
     # sets employing a clr transformation.  To obtain Mhalanobis distances 
     # a separate ilr transform of the data is undertaken with the degrees
     # of freedom being reduced by one.
     #
     # Note this procedure uses svd() rather than the classic solve().
     #
     # PCA output may be plotted with gx.rqpca.plot(), gx.rqpca.loadplot,
     # and gx.rqpca.screeplot(), and Mahalanobis distances may be plotted
     # with gx.md.plot().  The PCA solution may be rotated using gx.rotate().
     #
     if(!is.matrix(xx)) stop(deparse(substitute(xx)), " is not a Matrix")
     # Remove any rows containing NAs
     temp.x <- remove.na(xx)
     x <- temp.x$x; n <- temp.x$n; p <- temp.x$m
     # Save variable names and matrix row numbers
     matnames <- dimnames(xx)
     if(is.null(matnames[[1]])) matnames[[1]] <- c(1:n)
     wts <- numeric(n); wts[1:n] <- 1
     # Issue any warnings re small sample population size
     nc <- n
     cat("  n = ", n, "\tnc = ", nc, "\tp = ", p, "\t\tnc/p = ", round(nc/p, 2), "\n")
     if(nc < 5 * p)
         cat("  *** Proceed with Care, Core Size is < 5p ***\n")
     if(nc < 3 * p)
         cat("  *** Proceed With Great Care, nc = ", nc, ", which is < 3p ***\n")
     #
     # Perform a clr PCA
     x.clr <- clr(x, ifwarn = TRUE)
     dimnames(x.clr) <- matnames
     save <- cov.wt(x.clr, wt = wts, cor = TRUE)
     center <- save$center
     sd <- sqrt(diag(save$cov))
     corr <- save$cor
     temp <- sweep(x.clr, 2, center, "-")
     snd <- sweep(temp, 2, sd, "/")
     w <- sweep(temp, 2, sqrt(n) * sd, "/")
     wt <- t(as.matrix(w))
     a <- wt %*% as.matrix(w)
     b <- svd(a)
     cat("  Eigenvalues:", signif(b$d, 4), "\n")
     sumc <- sum(b$d)
     econtrib <- 100 * (b$d/sumc)
     cat("     as %ages:", round(econtrib, 1), "\n")
     rqscore <- w %*% b$v
     vcontrib <- numeric(p)
     for (j in 1:p) vcontrib[j] <- var(rqscore[, j])
     cat("  Score S^2s :", signif(vcontrib, 4), "\n")
     b1 <- b$v * 0
     diag(b1) <- sqrt(b$d)
     rload <- b$v %*% b1
     rcr <- rload[,  ] * 0
     rcr1 <- apply(rload^2, 1, sum)
     rcr <- 100 * sweep(rload^2, 1, rcr1, "/")
     dimnames(rload)[[1]] <- dimnames(rcr)[[1]] <- matnames[[2]]
     #
     # Compute Mahalanobis distances following an ilr transformation
     x.ilr <- ilr(x, ifwarn = FALSE)
     # Estimate ilr covariance matrix
     save.ilr <- cov.wt(x.ilr, wt = wts, cor = TRUE)
     # Invert ilr covariance matrix for use in gx.mvalloc.closed
     inverted <- ginv(save.ilr$cov)
     V <- orthonorm(p)
     cov.clr <- V %*% save.ilr$cov %*% t(V)
     dimnames(cov.clr)[[1]] <- dimnames(cov.clr)[[2]] <- matnames[[2]]
     inverted.clr <- V %*% inverted %*% t(V)
     dimnames(inverted.clr) <- dimnames(cov.clr)
     # Compute Mahalanobis distances and probabilities
     md <- mahalanobis(x.ilr, save.ilr$center, save.ilr$cov)
     p.ilr <- p - 1
     temp <- (nc - p.ilr)/(p.ilr * (nc + 1))
     ppm <- 1 - pf(temp * md, p.ilr, nc - p.ilr)
     epm <- 1 - pchisq(md, p.ilr)
     #
     invisible(list(main = main, input = deparse(substitute(xx)), proc = "cov",
         n = n, nc = nc, p = p, ifilr = TRUE, matnames = matnames, wts = wts,
         mean = center, cov = save$cov, cov.inv = inverted.clr, sd = sd,
         snd = snd, r = save$cor, eigenvalues = b$d, econtrib = econtrib,
         eigenvectors = b$v, rload = rload, rcr = rcr, rqscore = rqscore,
         md = md, ppm = ppm, epm = epm, nr = NULL))
}
