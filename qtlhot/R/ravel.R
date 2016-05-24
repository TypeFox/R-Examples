######################################################################
#
# ravel.R
#
# copyright (c) 20011, Brian S Yandell
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
# Part of the R/qtlhot package
# Contains: get.hotsize, get.hotscan,
#           ravel.scan,ravel.filter,get.hotravel
#
######################################################################
## Code from ~/p/private/diabetes1/diabetes10/Rcode/byandell/peaks/routines.R
##################################################################
get.hotsize <- function(x, quants)
{
  ## quants is preset for desired significance LOD threshold
  ##        and truncated at lowest acceptable value.
  ## x is set of LODs across all traits (say at a locus).
  
  ## Start from smallest LOD threshold and work back.
  count <- top <- 0
  lod <- min(quants)
  for(i in rev(seq(length(quants)))) {
    ## Find which x's are at or above quantile threshold.
    wh <- which(x >= quants[i])
    n.i <- length(wh)
    if(n.i) {
      if(n.i >= i) {
        ## Record if count is above index and break out.
        count <- n.i
        lod <- quants[i]
        top <- sum(x >= quants[1])
        break
      }
      ## Reduce to just those at least this large. Saves time.
      x <- x[wh]
    }
    else
      break
  }
  c(count, as.numeric(lod), top)
}
##################################################################
get.hotscan <- function(scan.out, quant, level = .95, lowest = 5)
{
  ## Slow, dumb way.
  
  ## Get LOD thresholds by hotspot size counts from 1 to lowest LOD allowed.
  quants <- quant[paste(100 * level, "%", sep = ""), ]
  quants <- quants[quants >= lowest]
  
  ## Loops are bad in R, but scan.out is > 1Gb.
  ## So, I take one row (= locus on genome) at a time.
  hotsize <- rep(NA, nrow(scan.out))
  hotlod <- rep(NA, nrow(scan.out))
  for(i in seq(nrow(scan.out))) {
    tmp <- get.hotsize(unlist(scan.out[i,-(1:2)]), quants)
    print(c(i,tmp))
    hotsize[i] <- tmp[1]
    hotlod[i] <- tmp[2]
  }
  hotscan <- scan.out[,1:4]
  hotscan[,3] <- hotsize
  hotscan[,4] <- hotlod
  names(hotscan)[3:4] <- c("count","LOD")

  hotscan
}
###########################################################
ravel.scan <- function(scan.out, threshold = 5)
{
  scans <- t(as.matrix(scan.out[,-(1:2)]))
  n.trait <- nrow(scans)
  wh <- which(scans >= threshold)
  out <- data.frame(trait = 1 + ((wh - 1) %% n.trait),
               pos = 1 + floor((wh - 1) / n.trait),
               lod = scans[wh])
  attr(out, "n.pos") <- ncol(scans)
  out
}
ravel.filter <- function(ravel, chr, loddrop = 1.5)
{
  n.pos <- attr(ravel, "n.pos")

  ## For each trait, keep only LODs within loddrop of peak
  ## on each chr.

  index <- paste(as.character(chr[ravel$pos]), ravel$trait, sep = ".")
  max.lod <- unlist(tapply(ravel$lod, index, max))
  ravel[(max.lod[index] - ravel$lod) <= loddrop, ]
}
###########################################################
get.hotravel <- function(ravel, quant, level = .95, lowest = 5,
                         verbose = FALSE)
{
  ## Better way?
  n.pos <- attr(ravel, "n.pos")
  
  ## Get LOD thresholds by hotspot size counts from 1 to lowest LOD allowed.
  quants <- quant[paste(100 * level, "%", sep = ""), ]
  quants <- quants[quants >= lowest]
  
  ## Loops are bad in R, but scan.out is > 1Gb.
  ## So, I take one row (= locus on genome) at a time.
  hottop <- hotlod <- hotsize <- rep(0, n.pos)
  for(i in seq(n.pos)) {
    tmp <- which(ravel$pos == i)
    if(length(tmp))
      tmp <- get.hotsize(ravel$lod[tmp], quants)
    if(verbose)
      print(c(i,tmp))
    hotsize[i] <- tmp[1]
    hotlod[i] <- tmp[2]
    hottop[i] <- tmp[3]
  }
  data.frame(count = hotsize, LOD = hotlod, singles = hottop)
}
