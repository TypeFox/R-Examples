######################################################################
# sim.R
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
# Contains: sim.null.cross, sim.null.pheno.data, include.hotspots,
#           mySimulations, sim.hotspot
######################################################################

## Generates a "null dataset" cross
##
sim.null.cross <- function(chr.len = rep(400,16), n.mar=185, n.ind = 112, type = "bc",
                           n.pheno = 6000, latent.eff = 1.5, res.var = 1,
                           init.seed = 92387475)
{
  set.seed(init.seed)
  mymap <- sim.map(len = chr.len, n.mar = n.mar, include.x=FALSE, eq.spacing=TRUE)
  scross <- sim.cross(map = mymap, n.ind = n.ind, type = type)
  scross <- calc.genoprob(scross, step=0)

  sim.null.pheno.data(cross = scross, n.pheno = n.pheno, latent.eff = latent.eff, res.var = res.var)
}
sim.null.pheno.data <- function(cross, n.pheno, latent.eff, res.var)
{
  n <- nind(cross)
  latent <- rnorm(n, 0, sqrt(res.var))
  ErrorM <- matrix(rnorm(n * n.pheno, 0, sqrt(res.var)), n, n.pheno)
  pheno <- data.frame(latent*latent.eff + ErrorM)
  names(pheno) <- paste("P", 1:n.pheno, sep="")
  cross$pheno <- pheno
  
  cross
}

## Generate hotspots for the simulated examples in the manuscript.
include.hotspots <- function(cross,
                             hchr,
                             hpos,
                             hsize,
                             Q.eff,
                             latent.eff,
                             lod.range.1,
                             lod.range.2,
                             lod.range.3,
                             res.var=1,
                             n.pheno,
                             init.seed)
{
  get.closest.pos.nms <- function(pos, cross, chr)
  {
    ## map <- attributes(cross$geno[[chr]]$prob)$map
    ## This can be simplified.
    map <- attributes(cross$geno[[chr]]$prob)$map
    ## map <- pull.map(cross, chr)
    q.nms <- names(map)
    map.pos <- as.numeric(map)
    tmp <- which.min(abs(map.pos-pos))
    closest.pos <- map.pos[tmp]
    nms <- q.nms[tmp]
    list(nms,closest.pos)
  }
  pull.prob <- function(cross)
  {
    out <- vector(mode = "list", length = nchr(cross))
    names(out) <- names(cross$geno)
    for(i in names(out))
      out[[i]] <- cross$geno[[i]]$prob  
    out
  }
  get.qtl.eff <- function(n, lod, res.var, latent.eff, Q.eff)
  {
##    lod <- runif(hsize, lod.range[1], lod.range[2])
##    r2 <- 1 - 10 ^ (-2 * lod / nind(cross))
##    beta <- sqrt(r2 * res.var * (1 + latent.eff ^ 2) / (Q.eff ^ 2 * (1 - r2) - res.var * r2))
    r2 <- 1 - 10^(-2*lod/n)
    sqrt(r2*(1 + latent.eff^2)/(Q.eff^2*(1 - r2) - r2))
  }
  update.pheno <- function(cross, hchr, hpos, hsize, Q.eff, latent.eff, lod.range,
                           res.var, index, hk.prob)
  {
    M.pos <- get.closest.pos.nms(hpos, cross, hchr)[[1]]
    M.dummy <- hk.prob[[hchr]][, M.pos, 1] - hk.prob[[hchr]][, M.pos, 2]
    M <- M.dummy * Q.eff + rnorm(nind(cross), 0, sqrt(res.var))

    ## QTL effect
    beta <- get.qtl.eff(n=nind(cross), lod=runif(hsize, lod.range[1], 
      lod.range[2]), res.var, latent.eff, Q.eff)

    for(j in seq(length(index)))
      cross$pheno[, index[j]] <- beta[j] * M + cross$pheno[, index[j]]

    cross
  }

  set.seed(init.seed)
  hk.prob <- pull.prob(cross)

  ## Why 50 for first, 500 for 2nd and 3rd?
  ## Why strange lod.range for 2nd?
  index1 <- sample(1:n.pheno, hsize[1], replace = FALSE)
  cross <- update.pheno(cross, hchr[1], hpos[1], hsize[1], Q.eff, latent.eff,
                        lod.range.1, res.var, index1, hk.prob)

  index2 <- sample((1:n.pheno)[-index1], hsize[2], replace = FALSE)
  cross <- update.pheno(cross, hchr[2], hpos[2], hsize[2], Q.eff, latent.eff,
                        lod.range.2, res.var, index2, hk.prob)

  index3 <- sample((1:n.pheno)[-c(index1, index2)], hsize[3], replace = FALSE)
  cross <- update.pheno(cross, hchr[3], hpos[3], hsize[3], Q.eff, latent.eff,
                        lod.range.3, res.var, index3, hk.prob)

  cross
}
#################################################################################
mySimulations <- function(...) sim.hotspot(...)
#################################################################################
sim.hotspot <- function(nSim, 
                        cross, 
                        n.pheno,
                        latent.eff,
                        res.var = 1,
                        n.quant,
                        n.perm,
                        alpha.levels,
                        lod.thrs,
                        drop.lod=1.5,
                        verbose = FALSE)
{
  s.quant <- seq(n.quant)

  nalpha <- length(alpha.levels)
  nlod <- length(lod.thrs)

  ## outputs count the number of times we detected
  ## a hotspot using the respective method
  outNL <- matrix(0, n.quant, nalpha)
  outN <- outWW <- matrix(0, nlod, nalpha)

  ## we are saving the thresholds of each simulation
  thrNL <- array(dim=c(n.quant, nalpha, nSim))
  thrN <- array(dim=c(nlod, nalpha, nSim))
  thrWW <- array(dim=c(nlod, nalpha, nSim))

  for(k in 1:nSim){
    mycat(k, verbose, TRUE)

    mycat("sim.null.pheno.data", verbose)
    ncross <- sim.null.pheno.data(cross, n.pheno, latent.eff, res.var)
  
    ## Simulate correlated phenotypes and create threshold summaries.
    out.sim <- filter.threshold(ncross, seq(n.pheno), latent.eff[k], res.var,
                             lod.thrs, drop.lod,
                             s.quant, n.perm, alpha.levels,
                             verbose)

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
