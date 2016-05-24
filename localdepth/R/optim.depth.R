#############################################################
#
#	optim.depth.simp.approx function
#	Author: Claudio Agostinelli and Mario Romanazzi
#	E-mail: claudio@unive.it
#	Date: August, 30, 2008
#	Version: 0.1-2
#
#	Copyright (C) 2008 Claudio Agostinelli and Mario Romanazzi
#
#############################################################

optim.depth.simp.approx <- function(x, y=NULL, tau, use=c('volume', 'diameter'), what=c('depth', 'localdepth'), weights=rep(1, nrow(y)), size=nrow(x)*4, size.final=nrow(x), subsize=0.1, nsamp=NULL, nsamp.final='all', q.level=0.4, tol=10^(-9), eps=10^(-2), verbose=FALSE, save=FALSE) {
  use <- match.arg(use)
  what <- match.arg(what)
  if (is.null(y))
    y <- x
  if (is.null(nsamp)) {
    nt <- choose(nrow(x), ncol(x)+1)
    if (nt < 500000) {
      nsamp <- 'all'
    } else {
      nsamp <- 500000
    }
  }
  ny <- nrow(y)
  y <- data.matrix(y)
  store <- list()
  store[[1]] <- NA
  maxdiff <- eps+1
  j <- 0
  while (maxdiff > eps) {
    j <- j + 1
    y <- rbind(sample.linearcomb(x=y, size=size, weights=weights, subsize=ceiling(subsize*nrow(y))), y)
#    if (j==1)
#      z <- rbind(z, y) 
    res <- localdepth.simp.approx(x=x, y=y, tau=tau, use=use, nsamp=nsamp, tol=tol)
    if (save) {
       store[[j]] <- res
    }
    if (what=='depth') {
      res <- res$depth
    } else {
      res <- res$localdepth
    }
#    if (j > 1) {
#      y <- rbind(z, y)
#      res <- c(res, weights)
#    } else
#      y <- z
##    qcut <- quantile(res, probs=q.level)
##    y <- y[res >=  qcut,]
    qcut <- ceiling(ny*q.level)
    resord <- order(res, decreasing=TRUE)
    y <- y[resord,][1:qcut,]
    
    maxdiff <- diff(range(res))/max(res)
    if (verbose) {
      cat('Max', max(res), '\n')
      cat('Diff', maxdiff, '\n')
      cat(summary(res), '\n')
    }
    weights <- res[resord][1:qcut]
  }
  ### final stage
  y <- rbind(sample.linearcomb(x=y, size=size.final, weights=weights, subsize=ceiling(subsize*nrow(y))), y)
  res <- localdepth.simp.approx(x=x, y=y, tau=tau, use=use, nsamp=nsamp.final, tol=tol)
  if (save) {
    res <- list(last=res, store=store)
  }
  return(res)
}
