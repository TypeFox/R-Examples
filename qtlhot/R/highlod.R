######################################################################
# highlod.R
#
# Elias Chaibub Neto
# Brian S Yandell
# Aimee Teo Broman
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
# Contains: highlod, sexbatch.covar, scanone.permutations,
#           cat.scanone, lod.quantile.permutation,
#           make.max.N, make.maxlod, smooth.neqtl
######################################################################
highlod <- function(scans, lod.thr = 0, drop.lod = 1.5,
                    extend = TRUE, restrict.lod = FALSE, ...)
{
  pheno.col <- seq(ncol(scans) - 2)

  if(is.null(lod.thr))
    lod.thr <- 0
  
  ## Extract matrix of lod scores.
  x <- as.matrix(scans[,-(1:2), drop = FALSE])

  ## Keep only traits with some LOD above lod threshold.
  keep <- apply(x, 2, function(x, lod.thr) any(x >= lod.thr), lod.thr)
  x <- x[, keep, drop = FALSE]
  
  ## Find which values are at within drop.lod of maximum per chr and trait.
  if(restrict.lod) {
    ## Restrict to loci above LOD threshold.
    if(extend)
      tmpfn <- function(x, lod.thr, drop.lod) {
        maxx <- max(x)
        g <- (maxx >= lod.thr) & (maxx <= x + drop.lod)
        if(any(g)) {
          d <- diff(g)
          ## Add one more pseudomarker on either side if possible.
          (g | (c(d,0) == 1) | (c(0,d) == -1)) & (x >= lod.thr)
        }
        else
          g
      }
    else
      tmpfn <- function(x, lod.thr, drop.lod) {
        (max(x) <= x + drop.lod) & (x >= lod.thr)
      }
  }
  else { ## Do not restrict support interval to be above lod.thr
    if(extend)
      tmpfn <- function(x, lod.thr, drop.lod) {
        maxx <- max(x)
        g <- (maxx >= lod.thr) & (maxx <= x + drop.lod)
        if(any(g)) {
          d <- diff(g)
          ## Add one more pseudomarker on either side if possible.
          g | (c(d,0) == 1) | (c(0,d) == -1)
        }
        else
          g
      }
    else
      tmpfn <- function(x, lod.thr, drop.lod) {
        maxx <- max(x)
        (maxx >= lod.thr) & (maxx <= x + drop.lod)
      }
  }
  lodint.pos <- function(x, chr, lod.thr, drop.lod) {
    unlist(tapply(x, chr, tmpfn, lod.thr, drop.lod))
  }
  wh <- apply(x, 2, lodint.pos, scans$chr, lod.thr, drop.lod)
  
  ## Get row and column indices.
  rr <- row(x)[wh]
  cc <- seq(keep)[keep][col(x)[wh]]

  ## Find which are within drop.lod of max lod per chr.
  lod <- x[wh]
  
  ## return data frame with genome row, trait column and lod value.
  out <- list(highlod = cbind.data.frame(row = rr, phenos = pheno.col[cc], lod = lod),
              chr.pos = scans[,1:2],
              names = names(scans)[-(1:2)])
  class(out) <- c("highlod", "list")
  attr(out, "lod.thr") <- lod.thr
  attr(out, "drop.lod") <- drop.lod
  out
}
print.highlod <- function(x, ...) print(summary(x, ...))
summary.highlod <- function(object, ...)
{
  summary(hotsize(object, ...))
}
plot.highlod <- function(x, ..., quant.level = NULL, sliding = FALSE)
{
  
  if(sliding) {
    if(is.list(quant.level))
      quant.level <- quant.level$max.lod.quant
    
    ## Need to supply quant.level as second argument.
    slidingbar.plot(slidingbar.create(x, quant.level, ...), ...)
  }
  else
    plot(hotsize(x, ..., quant.level = quant.level), ...)
}
###################################################################################
highlod.thr <- function(highobj, lod.thr)
{
  if(is.null(lod.thr))
      lod.thr <- attr(highobj, "lod.thr")

  if(!is.null(lod.thr)) {
    highobj$highlod <- highobj$highlod[highobj$highlod$lod >= lod.thr,, drop = FALSE]
    attr(highobj, "lod.thr") <- lod.thr
  }
  
  highobj
}
###################################################################################
cat.scanone <- function(dirpath = ".", filenames = permfiles, chr.pos)
{
  ## Folder should contain scanone highlods data across all traits for ONE permutation
  permfiles <- list.files(dirpath, paste("per.scan", "*", "RData", sep = "."))

  ## Make and remove per.scan.hl. Below use version from files.
  per.scan.hl <- NULL
  rm(per.scan.hl)
  
  for(i in 1:length(filenames)) {
    highobj <- with(filenames[i], {
      if(i==1) per.scan.hl else rbind.data.frame(highobj, per.scan.hl)
    })
  }
  cbind.data.frame(chr.pos[highobj$row,],highobj)
}
###################################################################################
sexbatch.covar <- function(cross, batch.effect, verbose = FALSE)
{
  ic <- getsex(cross)$sex
  ## Drop sex if only one present.
  if(length(unique(ic)) == 1)
    ic <- NULL
  
  if(!is.null(batch.effect)){
    batch <- cross$pheno[,batch.effect, drop = FALSE]
    tmp <- formula(paste("~ factor(", batch.effect, ")"))
    if(verbose)
      cat("sexbatch.covar", names(tmp), levels(factor(batch[[1]])), "\n")
    if(verbose)
      cat("sexbatch.covar", dim(batch), "\n")
    batch <- model.matrix(tmp,batch)[,-1, drop = FALSE]
    if(verbose)
      cat("sexbatch.covar", dim(batch), "\n")
    ac <- cbind(batch,ic)
  }
  else
    ac <- ic
  
  list(addcovar = ac, intcovar = ic)
}
## Performs and saves scanone on permuted dataset
scanone.permutations <- function(cross, pheno.col = seq(3, nphe(cross)),
                                 n.perm, seed=123456789, batch.effect = NULL,
                                 pheno.set = 1,
                                 lod.min, drop.lod = 1.5,
                                 addcovar = NULL, intcovar = NULL, ...)
{
  set.seed(seed[[1]])

  if(!is.null(batch.effect)) {
    cross <- subset(cross, ind = !is.na(cross$pheno[[batch.effect]]))
    covars <- sexbatch.covar(cross, batch.effect)
  }
  else
    covars <- list(addcovar = addcovar, intcovar = intcovar)
  
  n.ind <- nind(cross)

  perms <- matrix(NA, n.ind, n.perm)

  for(i in 1:n.perm){
    perm.cross <- cross
    perms[,i] <- tmp <- sample(c(1:n.ind), n.ind, replace=FALSE)
    perm.cross$pheno <- cross$pheno[tmp,]

    per.scan <- scanone(perm.cross, pheno.col=pheno.col, method="hk", 
                        addcovar=covars$addcovar, intcovar=covars$intcovar, ...)

    per.scan.hl <- highlod(per.scan, lod.thr = lod.min, drop.lod = drop.lod,
                                 restrict.lod = TRUE)$highlod

    save(per.scan.hl, perms,
         file=paste("per.scan",pheno.set, i,"RData",sep="."))
  }
}
###################################################################
pull.highlod <- function(object, chr, pos, ...)
{
  ## Kludge to get names if not in object
  pheno.names <- object$names
  if(is.null(pheno.names)) {
    extra <- list(...)
    m <- match("names", names(extra))
    if(!is.na(m))
      pheno.names <- extra[[m]]
  }
  wh.chr <- which(object$chr.pos$chr == chr)
  ## Want to expand this to handle range of positions...
  wh.pos <- which(object$chr.pos$pos[wh.chr] - min(pos) >= 0 &
                  object$chr.pos$pos[wh.chr] - max(pos) <= 0)
  wh.high <- which(object$highlod$row %in% wh.chr[wh.pos])
  wh.row <- object$highlod[wh.high, "row"]

  out <- data.frame(object$chr.pos[wh.row,], object$highlod[wh.high, -1])
  ## Add phenotype names if available.
  if(!is.null(pheno.names))
    out$phenos <- pheno.names[out$phenos]
  out
}
