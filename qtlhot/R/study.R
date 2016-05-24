######################################################################
# study.R
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
# Contains:  count.thr, exceed.thr,
#            get.hotspot, filter.threshold,
#            NL.counts, N.WW.counts, mycat
######################################################################

## This function computes the error rates for the NL-, N- and West/Wu methods
## out of nSim simulations. (It also returns the method's thresholds).
## At each simulation iteration this function: (1) generates a null dataset 
## cross; (2) perform Haley-Knott mapping for all traits; (3) applies the LOD
## drop interval computation to the scanone object; (4) determine the NL-, N-
## and West/Wu approches thresholds using the scanone results from step 3; 
## (5) for each of the three methods it computes the proportion of times out 
## of the nSim simulations we detected at least one false hotspot anywhere in
## the genome.  
##
count.thr <- function(scan, lod.thrs, droptwo = TRUE)
{
  ## Count number of traits at or above each value of lod.thrs for each locus.
  ## Result is n.loci * n.thr matrix.
  if(droptwo)
    scan <- scan[, -c(1,2), drop = FALSE]
  apply(scan, 1, exceed.thr, lod.thrs)
}
exceed.thr <- function(x, y)
{
  ## Finds out how many in x exceed each value of y.
  res <- rep(0, length(y))
  for(k in order(y)) {
    x <- x[x > y[k]]
    if(length(x) == 0)
      break
    res[k] <- length(x)
  }
  res
}
#################################################################################
get.hotspot <- function(filenames,
                        ## Following supplied in filenames[1].
                        n.quant, out.sim)
{
  ## filenames = list.files(".", paste(prefix, latent.eff, sets, "RData", sep = "."))
  ## latent.eff = 0, prefix = "Pilot", sets = "[0-9][0-9]*"
  
  ## Get stored n.quant value, and out.sim to get alpha.levels and lod.thrs.
  load(filenames[1])
  nSim <- length(filenames)
  n.quant <- n.quant ## Null action to make sure n.quant is valued.
  s.quant <- seq(n.quant)
  tmp <- dimnames(out.sim$N.thrs)
  lod.thrs <- as.numeric(tmp[[1]])
  alpha.levels <- as.numeric(tmp[[2]])

  ## May not have names on rows and columns.
  nalpha <- ncol(out.sim$N.thrs)
  nlod <- nrow(out.sim$N.thrs)
  
  ## outputs count the number of times we detected
  ## a hotspot using the respective method
  outNL <- matrix(0, n.quant, nalpha)
  dimnames(outNL) <- list(NULL, alpha.levels)
  outN <- outWW <- matrix(0, nlod, nalpha)
  dimnames(outN) <- dimnames(outWW) <- list(lod.thrs, alpha.levels)
  
  ## we are saving the thresholds of each simulation
  thrNL <- array(dim=c(n.quant, nalpha, nSim))
  thrN <- thrWW <- array(dim=c(nlod, nalpha, nSim))
  dimnames(thrN) <- dimnames(thrWW) <- list(lod.thrs, alpha.levels, NULL)
  dimnames(thrNL) <- list(NULL, alpha.levels, NULL)
  
  for(k in 1:nSim) {
    mycat(k, TRUE, TRUE)
    load(filenames[k])

    thrNL[,,k] <- out.sim$NL.thrs
    thrN[,,k] <- out.sim$N.thrs
    thrWW[,,k] <- out.sim$WW.thrs    
    outNL <- outNL + out.sim$NL
    outN <- outN + out.sim$N.counts
    outWW <- outWW + out.sim$WW.counts
  }
  
  NL.err <- outNL/nSim
  dimnames(NL.err) <- list(as.factor(s.quant), as.factor(alpha.levels))
  N.err <- outN / nSim
  dimnames(N.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  WW.err <- outWW / nSim
  dimnames(WW.err) <- list(as.factor(lod.thrs), as.factor(alpha.levels))
  list(nSim = nSim, NL.err=NL.err, N.err=N.err, WW.err=WW.err, thrNL=thrNL, thrN=thrN, 
       thrWW=thrWW)  
}
########################################################################################
filter.threshold <- function(cross, pheno.col, latent.eff, res.var,
                             lod.thrs, drop.lod = 1.5,
                             s.quant, n.perm, alpha.levels,
                             qh.thrs = summary(hotperm(cross, max(s.quant), n.perm, alpha.levels,
                               lod.thrs, verbose = verbose)),
                             ww.thrs = summary(ww.perm(highobj, n.perm, lod.thrs, alpha.levels)),
                             addcovar = NULL, intcovar = NULL,
                             verbose = FALSE, ...)
{
  mycat("scanone", verbose)
  scanmat <- scanone(cross, pheno.col = pheno.col, method = "hk", 
                     addcovar = addcovar, intcovar = intcovar, ...)

  ## Reduce to high LOD scores.
  mycat("highlod", verbose)
  highobj <- highlod(scanmat, min(lod.thrs), drop.lod, restrict.lod = TRUE)
  rm(scanmat)
  gc()
  
  ## computes an array of size n.quant by nalpha by npos.
  ## showing for each s.quant size and alpha level, the 
  ## hotspot sizes at each genomic location.
  mycat("NL.counts", verbose)
  NL.thrs <- qh.thrs[[1]]
  N.thrs <- qh.thrs[[2]]
  n.quant <- length(s.quant)
  NL <- NL.counts(highobj, n.quant, NL.thrs)
  
  ## computes a matrix of size nlod by npos.
  ## showing for each lod threshold the 
  ## hotspot sizes at each genomic location.
  mycat("N.WW.counts", verbose)
  N.WW <- N.WW.counts(highobj, lod.thrs, N.thrs, ww.thrs)
  
  list(NL.thrs = NL.thrs, N.thrs = N.thrs, WW.thrs = ww.thrs, NL = NL,
       N.counts = N.WW$N, WW.counts = N.WW$WW)
}

## Computes an array of size n.quant (number of spurious hotspots sizes) by 
## nalpha (number of significance levels) by npos (number of locus), and for
## each spurious hotspot size/significance level threshold, it computes the 
## number of traits mapping with LOD higher than the threshold at each one
## of the genomic positions.
##
NL.counts <- function(highobj, n.quant, NL.thrs)
{
  ## get the maximum spurious hotspot size (N-method) 
  ## for different QTL mapping significance levels

  XX <- quantile(highobj, n.quant = n.quant)
  NL.counts <- apply(NL.thrs, 2,
                     function(x,y) (x < y[seq(x)]),
                     XX)
  ## dimnames(NL.counts)[[2]] <- seq(n.quant)
  NL.counts
}

## Computes a matrix of size nlod (number of mapping thresholds) by npos 
## (number of locus), and for each LOD threshold, it computes the number
## of traits mapping with LOD higher than the threshold at each one of
## the genomic positions. The same counts are used by the N- and WW-methods.
##
N.WW.counts <- function(highobj, lod.thrs, N.thrs, WW.thrs)
{
  ## XX = genome position by number of traits above LOD threshold.
  XX <- max(highobj, lod.thr = lod.thrs)

  ## N.counts[lod,alpha] = TRUE if max hotspot size using lod is above alpha perm threshold.
  N.counts <- apply(N.thrs, 2, function(x,y) (x < y), XX)

  ## WW.counts[lod,alpha] = TRUE if max hotspot size using lod is above alpha perm threshold.
  WW.counts <- apply(WW.thrs, 2, function(x,y) (x < y), XX)
  dimnames(N.counts) <- dimnames(WW.counts) <- dimnames(N.thrs)

  list(N = N.counts, WW = WW.counts)
} 
mycat <- function(title, verbose = FALSE, init = FALSE, last = "\n")
{
  if(verbose) {
    if(verbose > 1) {
      if(init)
        cat("user system elapsed time\n")
      else
        cat(round(as.numeric(proc.time()[1:3])), "")
    }
    cat(title, last)
  }
}
