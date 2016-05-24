######################################################################
# westwu.R
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
# Contains: ww.perm, summary.ww.perm, print.ww.perm
######################################################################

## This function computes the West/Wu permutation thresholds.
## The output is a nlod (number of LOD thresholds) by nalpha (number of 
## significance levels) matrix, where each entry shows the hotspot size 
## significance threshold of the West/Wu approach. Note we have two "alphas" 
## here, one for the QTL mapping (the LOD thresholds) and one for the 
## permutation significance (alpha levels of lod.thrs).
##
## Note that I separated the original ww.permutations() into a piece that do the 
## actual permutations [ww.perm.matrix() function] and a piece that summarizes it
## [the ww.summary() function] in the same way you did with the NL.N.permutations()
## function.  
## 
ww.perm <- function(highobj, n.perm, lod.thrs, alpha.levels, verbose = FALSE)
  ww.perm.highlod(highobj, n.perm, lod.thrs, alpha.levels, verbose)
####################################################################################
print.ww.perm <- function(x, ...) print(summary(x, ...))
summary.ww.perm <- function(object, alpha.levels = attr(object, "alpha.levels"), ...)
{
  nalpha <- length(alpha.levels)
  ww.thrs <- t(apply(object, 2, quantile, 1 - alpha.levels))
  dimnames(ww.thrs) <- list(dimnames(object)[[2]], as.character(alpha.levels))
  ww.thrs
}
####################################################################################
ww.perm.highlod <- function(highobj, n.perm, lod.thrs, alpha.levels, verbose = FALSE)
{
  ## These permutations differ from Elias's original (see inst/deprecated.R):
  ## 1. Only sampling on subset of phenotypes with peaks.
  ##    Thus the number of random sample calls is smaller.
  ## 2. chr.pos -> sample(chr.pos) rather than lod <- sample(lod).
  ##    Thus shuffle in inverted: s <- sample(1:10) vs. o <- order(s)
  ## 3. Original code shuffled scanmat repeatedly
  ##    rather than new shuffle each time. This is a minor bug rather than a feature.
  ## 4. New code is much faster.

  ## Both versions break up the support intervals with shuffling.
  ## This can lead to more peaks across genome.
  ## ww.perm.maxlod (below) resolves this by focusing only on peaks,
  ## but only makes sense if ww.perm is used with window = 0.
  
  highobj <- highlod.thr(highobj, min(lod.thrs))
  
  phenos <- sort(unique(highobj$highlod$pheno))
  n.chrpos <- nrow(highobj$chr.pos)

  nlod <- length(lod.thrs)
  max.ww <- matrix(NA, n.perm, nlod)
  dimnames(max.ww) <- list(NULL, as.character(lod.thrs))
  
  ## Want to replace row in highobj$highlod with permuted row.
  ## With separate permutation by pheno of seq(n.chrpos).

  ## Shuffle rows of chr.pos. Done separately for each pheno using tapply below.
  mysam <- function(row, n) sample(n)[row]

  ## Order highlod by pheno to make it easier
  highobj$highlod <- highobj$highlod[order(highobj$highlod$pheno), ]
  highs <- highobj
  
  for(i in 1:n.perm) {
    if(verbose)
      cat("\n", i, "")

    ## Permute columns of the scan object separately by trait.
    mycat("sample...", verbose, last = "")
    row.perm <- unlist(tapply(highobj$highlod$row, highobj$highlod$pheno, mysam, n.chrpos))
    highs$highlod$row <- row.perm
    
    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("max...", verbose, last = "")
    tmp <- max(highs, lod.thr = lod.thrs)
    if(length(tmp)) {
      m <- match(tmp$lod.thr, lod.thrs)
      max.ww[i, m] <- tmp$max.N
    }
  }
  if(verbose) cat("\n")
  
  class(max.ww) <- c("ww.perm", "list")
  attr(max.ww, "lod.thrs") <- lod.thrs
  attr(max.ww, "alpha.levels") <- alpha.levels
  
  max.ww
}
####################################################################################
ww.perm.maxlod <- function(highobj, n.perm, lod.thrs, alpha.levels, verbose = FALSE)
{
  highobj <- highlod.thr(highobj, min(lod.thrs))
  maxobj <- pull.max(highobj)
  
  phenos <- sort(unique(highobj$highlod$pheno))
  n.chrpos <- nrow(highobj$chr.pos)

  nlod <- length(lod.thrs)
  max.ww <- matrix(NA, n.perm, nlod)
  dimnames(max.ww) <- list(NULL, as.character(lod.thrs))
 
  ## Shuffle rows of chr.pos. Done separately for each pheno using tapply below.
  mysam <- function(row, n) sample(n)[row]
 
  ## Order highlod by pheno to make it easier
  maxobj$highlod <- maxobj$highlod[order(maxobj$highlod$pheno), ]
  highs <- maxobj
  
  for(i in 1:n.perm) {
    if(verbose)
      cat("\n", i)

    ## Permute columns of the scan object separately by trait.
    mycat("sample...", verbose, last = "")
    row.perm <- unlist(tapply(maxobj$highlod$row, maxobj$highlod$pheno, mysam, n.chrpos))
    highs$highlod$row <- row.perm
    
    ## Get the maximum spurious hotspot size (N-method) across genome
    ## for different QTL mapping significance levels.
    mycat("max...", verbose, last = "")
    tmp <- max(highs, lod.thr = lod.thrs)
    if(length(tmp)) {
      m <- match(tmp$lod.thr, lod.thrs)
      max.ww[i, m] <- tmp$max.N
    }
  }
  class(max.ww) <- c("ww.perm", class(max.ww))
  attr(max.ww, "lod.thrs") <- lod.thrs
  attr(max.ww, "alpha.levels") <- alpha.levels
  
  max.ww
}
######################################################################
pull.max <- function(highobj)
{
  highobj$highlod <- highobj$highlod[order(highobj$highlod$pheno), ]
  chr.pheno <- ordered(paste(highobj$chr.pos$chr[highobj$highlod$row],
                             highobj$highlod$pheno, sep = "."))
  tmpfn <- function(x) sample(which.max(x), 1)
  tmp <- unlist(tapply(highobj$highlod$lod, chr.pheno, which.max))
  tmp2 <- table(chr.pheno)
  highobj$highlod <- highobj$highlod[tmp + c(0, cumsum(tmp2)[-length(tmp2)]), ]
  highobj
}




