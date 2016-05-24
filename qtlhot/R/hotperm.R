######################################################################
# perm.R
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
# Contains: hotperm,summary.hotperm,print.hotperm
######################################################################

## This set of functions compute the permutation LOD thresholds for the NL-method
## and the permutation hotspot size thresholds for the N-method. The output
## is a list with two elements: the NL- and the N-method's threshold matrices.
## The NL-method output is a nN (number of spurious hotspot sizes) by nalpha
## (number of significance levels) threshold matrix. Note that for the
## NL-method we have a single "alpha" since we use the same significance level
## for QTL mapping and permutation significance.
## The N-method output is a nlod (number of LOD thresholds) by nalpha (number
## of significance levels) threshold matrix. Note that here we have two
## "alphas", one for the QTL mapping (the LOD thresholds) and one for the 
## permutation significance (alpha levels). 
##
hotperm <- function(cross, n.quant, n.perm, lod.thrs, alpha.levels, drop.lod = 1.5,
                    window = NULL, verbose = FALSE, init.seed = 0,
                    addcovar = NULL, intcovar = NULL, ...) 
{
  set.seed(init.seed)
  n.phe <- nphe(cross)
  pheno.col <- seq(n.phe)
  n.ind <- nind(cross)
  n.quant <- min(n.quant, n.phe)

  tmp <- table(sapply(cross$pheno, class))
  if(length(tmp) > 1 | names(tmp)[1] != "numeric")
    stop("all phenotypes in cross object must be numeric")
  
  s.quant <- seq(n.quant)
  quants <- 1 - (s.quant - 1) / n.phe
  n.lod <- length(lod.thrs)

  max.N <- matrix(0, n.perm, n.lod)
  dimnames(max.N) <- list(NULL, as.character(lod.thrs))
  if(is.null(window))
    max.N.window <- NULL
  else
    max.N.window <- max.N

  max.lod.quant <- matrix(0, n.perm, n.quant)
  dimnames(max.lod.quant) <- list(NULL, as.character(s.quant))

  for(i in 1:n.perm){
    if(verbose)
      cat("\n", i, "")
    
    ## permute rows of the phenotype data matrix
    perm.cross <- cross
    tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    ## perform mapping analysis in the permuted data
    mycat("scanone...", verbose, last = "")
    ## NB: scanone groups phenos in batches based on missing data patterns.
    scanmat <- scanone(perm.cross, pheno.col = pheno.col, method = "hk", 
                        addcovar = addcovar, intcovar = intcovar, ...)

    ## Reduce to high LOD scores.
    mycat("highlod...", verbose, last = "")
    highs <- highlod(scanmat, min(lod.thrs), drop.lod, restrict.lod = TRUE)
    rm(scanmat)
    gc()

    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("max...", verbose, last = "")
    maxhi <- max(highs, lod.thr = lod.thrs, window = window)
    max.N[i, ] <- maxhi$max.N
    if(!is.null(window))
      max.N.window[i, ] <- maxhi$max.N.window
    
    ## get the maximum lod-quantile across the genome
    ## rows indexes the permutations
    ## columns indexes the s.quant quantiles
    mycat("quantile...", verbose, last = "")
    tmp <- quantile(highs, n.quant = n.quant)
    if(length(tmp))
      max.lod.quant[i, seq(tmp)] <- tmp
  }
  if(verbose) cat("\n")
  
  out <- list(max.N = max.N, max.N.window = max.N.window,
              max.lod.quant = max.lod.quant)
  class(out) <- c("hotperm", "list")
  attr(out, "lod.thrs") <- lod.thrs
  attr(out, "alpha.levels") <- alpha.levels
  
  out
}
print.hotperm <- function(x, ...) print(summary(x, ...))
summary.hotperm <- function(object, quant.levels, ...)
{
  out <- quantile(object, ...)

  attr(out, "window") <- attr(object, "window")

  alpha.levels <- attr(object, "alpha.levels")
  if(max(alpha.levels) < 0.5)
    alpha.levels <- 1 - alpha.levels
  n.quant <- ncol(object$max.lod.quant)
  if(missing(quant.levels)) {
    quant.levels <- log10(n.quant)
    quant.levels <- round(10 ^ c(outer(log10(c(1,2,5)), seq(0, floor(quant.levels)), "+")))
  }
  quant.levels <- quant.levels[quant.levels <= n.quant]
  if(max(quant.levels) < n.quant)
    quant.levels <- c(quant.levels, n.quant)
  out$max.lod.quant <- t(apply(object$max.lod.quant[, quant.levels, drop = FALSE],
                               2, quantile, probs = alpha.levels, na.rm = TRUE))
  
  class(out) <- c("summary.hotperm", class(out))
  out
}
print.summary.hotperm <- function(x, ...)
{
  cat("max.N: hotspot threshold by single-trait LOD threshold and significance level\n")
  print(ceiling(x$max.N))
  window <- attr(x, "window")
  if(!is.null(window)) {
    cat(paste("\nmax.N.window: smoothed hotspot threshold by single-trait LOD threshold and significance level ",
              "(window = ", window, ")\n", sep = ""))
    print(ceiling(x$max.N.window))
  }
  cat("\nmax.lod.quant: LOD threshold by hotspot size quantile and significance level\n")
  print(round(x$max.lod.quant, 2))
  invisible()
}
plot.hotperm <- function(x, probs = seq(0.9, 0.99, by = 0.01), level = 0.95, ...)
{
  lod.thrs <- attr(x, "lod.thrs")
  alpha.levels <- attr(x, "alpha.levels")
  if(max(alpha.levels) <= 0.5)
    alpha.levels <- 1 - alpha.levels
  if(max(probs) <= 0.5)
    probs <- 1 - probs
  if(level < 0.5)
    level <- 1 - level
  
  lod.thr <- lod.thrs[which.min(abs(level - alpha.levels))[1]]
  out <- quantile(x, probs, lod.thr = lod.thr, ...)

  wh.thr <- which.min(abs(lod.thr - lod.thrs))[1]
  wh.level <- which.min(abs(level - probs))[1]
  
  tmp.plot <- function(x.vals, quant, x.crit, probs, level, wh.thr, is.quantile = FALSE, main = "",
                       add.level = FALSE)
  {
    n.probs <- length(probs)
    wh.level <- which.min(abs(level - probs))[1]
    quant.thr <- quant[wh.thr, wh.level]

    xlabs <- "single trait LOD threshold"
    if(is.quantile)
      xlabs <- paste(xlabs, "quantile")
    
    plot(range(x.vals), c(0, max(quant)), type = "n", xlab = "", ylab = "")
    mtext(xlabs, 1, 2)
    mtext("hotspot size", 2, 2)
    abline(v = x.crit, col = "darkgray", lty = 2)
    abline(h = quant.thr, col = "darkgray", lty = 2)
    mtext(ceiling(quant.thr), 2, at = quant.thr, las = 2, cex = 1)
    for(i in seq(n.probs)) {
      lines(rev(sort(x.vals)), quant[,i],
            lwd = 1 + 2 * (round(probs[i] - level, 2) == 0))
    }

    text(x.crit, quant[wh.thr, n.probs] + 5, 1 - max(probs), adj = 0)
    text(x.crit, quant[wh.thr, 1] - 5, 1 - min(probs), adj = 1)
    text(par("usr")[1], quant.thr + 5, 1 - level, adj = 0)

    if(add.level)
      main <- paste(main, "\n hotspot size significance level =", 1 - max(probs), "to", 1 - min(probs))
    mtext(main, 3, 0.5)
  }

  tmpar <- par(mfrow = c(1 + !is.null(out$max.N.window),2), mar = c(3.1,3.1,2.1,0.1))
  if(!is.null(out$max.N.window)) {
    ## Jansen method, smoothing.
    tmp.plot(lod.thrs, out$max.N.window, lod.thr, probs, level, wh.thr, FALSE,
             "Jansen method 5cM window")
    
    tmp.plot(probs, out$max.N.window, level, probs, level, wh.thr, TRUE,
             "Jansen method 5cM window")
  }

  tmp.plot(lod.thrs, out$max.N, lod.thr, probs, level, wh.thr, FALSE,
           "Jansen method per locus")
  tmp.plot(probs, out$max.N, level, probs, level, wh.thr, TRUE,
           "Jansen method per locus")
  par(tmpar)

  if(!is.null(out$max.lod.quant)) {
    n.quant <- nrow(out$max.lod.quant)
    n.probs <- length(probs)    
    quant.thr <- max(which(out$max.lod.quant[,wh.level] >= lod.thr), na.rm = TRUE)
    
    ## Chaibub Neto method.
    plot(c(1,n.quant), range(out$max.lod.quant, na.rm = TRUE), type = "n",
         xlab = "significant hotspot size with given threshold",
         ylab = "hotspot LOD score threshold",
         log = "xy")
    abline(h = lod.thr, col = "darkgray", lty = 2)
    abline(v = quant.thr, col = "darkgray", lty = 2)
    mtext(ceiling(quant.thr), 1, at = quant.thr, las = 2)
    for(i in seq(n.probs)) {
      if(any(!is.na(out$max.lod.quant[,i])))
        lines(seq(n.quant), out$max.lod.quant[,i],
              lwd = 1 + 2 * (round(probs[i] - level, 2) == 0))
    }
    
    n.thr2 <- length(lod.thrs) / 2
    text(n.thr2 + 1, out$max.lod.quant[n.thr2, n.probs], 1 - max(probs), adj = 0)
    text(n.thr2 - 1, out$max.lod.quant[n.thr2, 1], 1 - min(probs), adj = 1)
    text(par("usr")[1], lod.thr, level, adj = 0)
    title(paste("hotspot LOD threshold by hotspot size\nsignificance level =",
                1 - max(probs), "to", 1 - min(probs)))
  }
}
