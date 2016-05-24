######################################################################
# diag.R
#
# Brian S Yandell
#
#     This program is free software; you can redistribute it and/or
#     modify it under the terms of the GNU General Public License,
#     version 3, as published by the Free Software Foundation.
# 
#     This program is distributed in the hope that it will be useful,
#     but without any warranty; without even the implied warranty of
#     merchantability or fitness for a particular purpose.  See the GNU
#     General Public License, version 3, for more details.
# 
#     A copy of the GNU General Public License, version 3, is available
#     at http://www.r-project.org/Licenses/GPL-3
#
# Contains: dist.qtlnet, edgematch.qtlnet, mds.qtlnet, plotbic.qtlnet,
# 
######################################################################

dist.qtlnet <- function(qtlnet.object, min.prob = 0.9, method = "manhattan", cex = 5)
{
  ## Fold to unique edges; threshold on min.prob.
  M <- apply(qtlnet.object$Mav, 3, function(x) 1 * (fold.M(x) >= min.prob))
  
  mbic <- meanbic(qtlnet.object)
  wh <- which.min(mbic)
  
  out <- list(sum = apply(M, 2, sum),
              common = apply(M, 2, function(x,y) sum(x*y), M[,wh]),
              wh = wh)
  
  plot(jitter(out$sum), jitter(out$common), xlab = "edges", ylab = "common with best")
  abline(0,1,col = "gray")
  abline(h = out$sum[wh], v = out$sum[wh], col = "gray")
  out
}
######################################################################
edgematch.qtlnet <- function(qtlnet.object, min.prob = 0.9, method = "manhattan", cex = 5)
{
  ## Deconvolve this fold to find out what pairs.
  ## Fold to unique edges; threshold on min.prob.
  M <- apply(qtlnet.object$Mav, 3, function(x) 1 * (abs(fold.M(x)) >= min.prob))
  
  mbic <- meanbic(qtlnet.object)
  wh <- which.min(mbic)
  common <- apply(M, 2, function(x,y) x*y == 1, M[,wh])
  extra <- apply(M, 2, function(x,y) x*(1-y) == 1, M[,wh])
  
  plot(c(1, nrow(common)), c(1, ncol(common)), type = "n", xlab = "edge", ylab = "run")
  abline(v = which(common[,wh]), col = "black")
  for(i in seq(ncol(common))) {
    if(i != wh & any(extra[,i]))
      points(which(extra[,i]), rep(i, sum(extra[,i])), col = "gray")
    points(which(common[,i]), rep(i, sum(common[,i])), col = ifelse(i==wh, "red", "black"))
  }
}
######################################################################
mds.qtlnet <- function(qtlnet.object, min.prob = 0.9, method = "manhattan", cex = 5)
{
  M <- apply(qtlnet.object$Mav, 3, fold.M)
  d <- dist(t(M), method) # euclidean distances between the cols
  fit <- cmdscale(d,eig=TRUE, k=2) # k is the number of dim
  
  ## plot solution
  x <- fit$points[,1]
  y <- fit$points[,2]
  mbic <- meanbic(qtlnet.object)
  wh <- which.min(mbic)
  rbic <- range(mbic)
  cex <- 1 + (cex - 1) * (rbic[2] - mbic) / diff(rbic)
  plot(x, y, xlab="Coordinate 1", ylab="Coordinate 2",
       main="Metric MDS", type="p", cex = cex)
  points(x[wh], y[wh], cex = cex[wh], col = "red")
  ## text(x, y, labels = row.names(M), cex=.7)
  invisible(list(M=M, d=d, fit = fit, mbic = mbic))
}
######################################################################
plotbic.qtlnet <- function(x, ..., smooth = TRUE)
{
  nSamples <- attr(x, "nSamples")
  runs <- length(nSamples)
  burnin <- attr(x, "burnin")
  
  rngfn <- function(x, burnin) {
    tmp <- which(seq(x) >= burnin * length(x))
    range(x[tmp])
  }
  plotfn <- function(post.bic, burnin, col = "gray") {
    tmp <- which(seq(post.bic) >= burnin * length(post.bic))
    lines(tmp, post.bic[tmp], col = col)
  }
  splotfn <- function(post.bic, burnin) {
    tmp <- which(seq(post.bic) >= burnin * length(post.bic))
    lines(tmp, lowess(post.bic[tmp])$y, col = "black")
  }

  bicol <- ifelse(smooth, "gray", "black")
  if(runs == 1) {
    range.bic <- rngfn(x$post.bic, burnin)
    plot(c(1,max(nSamples)),range.bic, type = "n",
         xlab = "Sample Index", ylab = "BIC")
    plotfn(x$post.bic, burnin, bicol)
  }
  else {
    run.id <- rep(seq(runs), nSamples)
    
    range.bic <- range(unlist(tapply(x$post.bic, run.id,
                                     rngfn, burnin)))
    plot(c(1,max(nSamples)),range.bic, type = "n",
         xlab = "Sample Index", ylab = "BIC")
    
    tapply(x$post.bic, run.id, plotfn, burnin, bicol)
    if(smooth)
      tapply(x$post.bic, run.id, splotfn, burnin)
  }
  title(paste("BIC samples for", runs, "MCMC", ifelse(runs == 1, "run", "runs")))
}

######################################################################
newfun <- function(qtlnet.object, burnin = attr(qtlnet.object, "burnin"),
                   wh = which.min(meanbic(qtlnet.object, burnin)))
{
  ## Nice idea, but not working the way I thought.
  ## Want sumM to be score of posterior for edge.
  M1 <- qtlnet.object$M[,,wh]
  M1 <- t(apply(M1, 1,
                function(x) {
                  s <- sum(x)
                  if(s > 0)
                    x <- x / s
                  x
                }))
  sumM <- M2 <- M1
  for(i in seq(2, nrow(M1) - 1)) {
    M2 <- M1 %*% M2
    sumM <- sumM + M2
  }
  upM <- (sumM + t(sumM))[upper.tri(sumM)]

  runM <- apply(qtlnet.object$M, 1:2, sum)
  runM <- (runM + t(runM))[upper.tri(runM)]
  plot(upM,runM)
  M <- qtlnet.object$M[,,wh]
  whs <- which(M[upper.tri(M)] > 0.9 | t(M)[upper.tri(M)] > 0.9)
  points(upM[whs], runM[whs], col = "red")

  ## This is easier to understand.
  ## Does M1 have an edge, and does it agree with most of the runs?
  M1u <- (M1+t(M1))[upper.tri(M1)]
  plot(M1u, runM)
  points(M1u[whs], runM[whs], col = "red")
  abline(v=0.9)
}
zero.M <- function(qtlnet.object, run = which.min(mbic),
                   burnin = attr(qtlnet.object, "burnin"))
{
  ## apply(out.qtlnet$Mav,3, function(x) sum(x >.9))
  ## round(apply(out.qtlnet$Mav,3, function(x) mean(x[x >.9])),3)
  nSamples <- attr(qtlnet.object, "nSamples")
  runs <- length(nSamples)
  run.id <- rep(seq(runs), nSamples)
  M0 <- attr(qtlnet.object, "M0")
  ravel <- row(as.matrix(M0)) > col(as.matrix(M0))
  
  tmpfn <- function(post.model, burnin, ravel) {
    M <- apply(model2M(post.model), 1:2, sum)
    (M[ravel] + t(M)[ravel]) == 0
  }
  cat("Extracting network matrices...\n")
  out <- matrix(unlist(tapply(qtlnet.object$post.model, run.id, tmpfn,
                              burnin, ravel)),
                ncol = runs)
  
  mbic <- meanbic(qtlnet.object, burnin)
  wh <- which.min(mbic)
  
  data.frame(nonzero = apply(out, 2, sum),
             agree = apply(out, 2, function(x,y) sum(x == y & y > 0),
               out[,run]),
             mean.bic = mbic)
}
best.qtlnet <- function(x, burnin = attr(x, "burnin"),
                        wh = which.min(meanbic(x, burnin)))
{
  subset(x, wh)
}
meanbic <- function(qtlnet.object, burnin = attr(qtlnet.object, "burnin"))
{
  nSamples <- attr(qtlnet.object, "nSamples")
  runs <- length(nSamples)
  run.id <- rep(seq(runs), nSamples)
  
  tmpfn <- function(x, burnin) {
    tmp <- which(seq(x) >= burnin * length(x))
    mean(x[tmp])
  }
  mbic <- tapply(qtlnet.object$post.bic, run.id, tmpfn, burnin)
}
######################################################################
fold.M <- function(x) {
  low <- upper.tri(x)
  t(x)[low] - x[low]
}
######################################################################


  
