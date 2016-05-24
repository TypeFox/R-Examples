#
# File: SGSE.R
# Author/Maintainer: rob.frost@dartmouth.edu
# Description: Depends on the PCGSE.R for the pcgse() function. Implements spectral gene set enrichment (SGSE) method, which
#              computes the statistical enrichment of one or more variable groups on the spectra of a empirical data set. 
#              Spectral enrichment is calculated by combining the PCGSE p-values for each gene set using the weighted Z-method 
#              weights equal to the lower-tailed p-values associated with the Gamma-approximation of the
#              Tracey-Widom Law of Order 1 distributed centered and scaled eigenvalues. 
#                
# Copyright (C) Dartmouth College
#

library(RMTstat)

#-------------------------------------------------------------------------------------------------------------------------------
# Public methods
#-------------------------------------------------------------------------------------------------------------------------------

#
# Compute the statistical enrichment or depletion of gene sets on the significant PC of the specified data set. 
# The enrichment on each PC is computed using the pcgse() method. The number of PCs to test is determined by the
# pc.selection.method parameter. A single overall enrichment p-value is computed by combining the individual PC 
# enrichment p-values using Stouffer's method with the -log(PC p-value) from the RMT-based test as the weight.
#
# Inputs:
#   All of the input parameters accepted by the pcgse() function in the PCGSE package:
#      data, prcomp.output, gene.sets, gene.statistic, transformation, gene.set.statistic,
#      gene.set.test, nperm
#   pc.selection.method: Method used to determine the PCs for which enrichment will be computed. One of the following:
#         "all": All PCs with non-zero variance will be used.
#         "specific": The set of PCs specified by pc.indexes will be used.
#         "rmt": The set of PCs with significant eigenvalues according to the Tracy-Widom distribution for a white Wishart at the specified alpha.
#   pc.indexes: Indices of the PCs for which enrichment should be computed. Must be specified if pc.selection.method is "specific".
#   rmt.alpha: Significance level for selection of PCs according to the Tracy-Widom distribution. Must be specified if pc.selection.method is "rmt".
#   pcgse.weight: Weight used with the weighted Z-method to combine the p-values from the PCGSE tests on all selected PCs for a specific gene set. 
#         Must be one of the following:
#         "variance": The PC variance is used as the weight. NOTE: this should only be used for evaluation and testing.
#         "rmt.scaled.var": Product of the PC variance and the Tracey-Widom lower-tailed p-value for the eigenvalue associated with the PC is used as the weight. 
# Output: 
#   List with the following elements:
#      pc.indexes: Indices of the PCs on which enrichment was performed.
#      pcgse: Output from pcgse() on the PCs identified by pc.indexes.
#      sgse: Vector of combined p-values for all PCs identified by pc.indexes.
#      weights: Vector of PC-specific weights for the PCs in pc.indexes.
#
sgse = function(data, prcomp.output=NA, gene.sets, 
    gene.statistic="z", transformation="none", 
    gene.set.statistic="mean.diff", gene.set.test="cor.adj.parametric", nperm=999,
    pc.selection.method="all", pc.indexes=NA, rmt.alpha=.05, pcgse.weight="rmt.scaled.var") {

  current.warn = getOption("warn")
  options(warn=-1)
  if (!(pc.selection.method %in% c("specific", "all", "rmt"))) {
    stop("pc.selection.method must be one of 'specific', 'rmt' or 'all'")
  }  
  if (!(pcgse.weight %in% c("variance", "rmt.scaled.var"))) {
    stop("pcgse.weight must be one of 'rmt.scaled.var' or 'variance'")
  }    
  if (is.na(data)) {
    stop("'data must' be specified!")
  }
  if (is.na(gene.sets)) {
    stop("gene.sets must be specified!")
  }
  # Compute PCA if necessary 
  if (is.na(prcomp.output)) {    
    prcomp.output=prcomp(data, center=T, scale=T)
  }      
  
  num.gene.sets = length(gene.sets)
  if (is.matrix(gene.sets)) {
    num.gene.sets = nrow(gene.sets)
  }
    
  # If needed, get the Tracey-Widom statistics for the PCs and the equivalent p-values
  pc.rmt = NA
  if (pc.selection.method == "rmt" | pcgse.weight == "rmt.scaled.var") {
    pc.rmt = pcRMT(prcomp.output)
  }
    
  # Determine which PCs to analyze  
  if (pc.selection.method == "rmt") {
    # analyze PCs whose RMT-based p-value is less than the specified alpha
    pc.indexes = which(pc.rmt$pvals < rmt.alpha)  
  } else if (pc.selection.method == "specific") {
    # analyze the specified PCs
    pc.indexes = pc.indexes
  } else if (pc.selection.method == "all") {
    # analyze all PCs with non-zero eigenvalues
    pc.indexes = which(prcomp.output$sdev > 0)
  }
  
  # Compute enrichment for the desired PCs
  pce = pcgse(data=data, prcomp.output=prcomp.output, pc.indexes=pc.indexes,
      gene.sets=gene.sets, gene.statistic=gene.statistic, transformation=transformation,
      gene.set.statistic=gene.set.statistic, gene.set.test=gene.set.test, nperm=nperm)
  
  # Compute combined PC enrichment using the weighted Z-method using either the 
  # Tracey-Widom scaled variance of the variance as the weight
  if (pcgse.weight == "variance") {
    weights=prcomp.output$sdev[pc.indexes]^2    
  } else if (pcgse.weight == "rmt.scaled.var") {
    weights=prcomp.output$sdev[pc.indexes]^2 * (1-pc.rmt$pvals[pc.indexes])  
  }
  sgse.pvalues = combinePValues(p.values=pce$p.values, method="weighted.Z",weights=weights)
  
  results = list(pc.indexes=pc.indexes, pcgse=pce, sgse=sgse.pvalues, weights=weights)  
  options(warn=current.warn) 
  
  return (results)
}

#-------------------------------------------------------------------------------------------------------------------------------
# Internal methods - RMT and p-value combination
#-------------------------------------------------------------------------------------------------------------------------------

#
# Center and scale the PC eigenvalues so that they would have be Tracey-Widom Law of Order 1 distribution under the
# null of a white Wishart matrix. The p-value can be computes using either the TW distribution from RMTstat or
# the Gamma approximation as outlined in http://arxiv.org/pdf/1209.3394v3.pdf.
#
# Inputs:
#  prcomp.output: Results from prcomp()
#
# Output:
#  List with two elements:
#   pvals: p-values
#   tws: median-centered TW statistics
#
pcRMT = function(prcomp.output) {
  eigenvals = prcomp.output$sdev^2
  r = nrow(prcomp.output$rotation)
  n = nrow(prcomp.output$x)
  pvals = rep(1, length(eigenvals))
  tws = rep(1, length(eigenvals))  
  for (i in 1:length(eigenvals)) {
    eval = eigenvals[i]
    pdim = (r-i+1)
    #tw.median = qtw(.5)#qWishartMax(.5, ndf=n, pdim=(r-i+1))
    par = WishartMaxPar(ndf=n, pdim=pdim)
    tws[i] = (eval-par$centering)/par$scaling 
    #pvals[i] = ptw(tws[i], lower.tail=F)
    # Not using the RMTstat p-value for TW RV since this is based on a rather sparse look-up table and
    # generates 0 for many large eigenvalues. Instead, use the gamma approximation (values of the
    # offset, shape and scale parameters are taken from Table 1 in http://arxiv.org/pdf/1209.3394v3.pdf. 
    pvals[i] = pgamma(tws[i]+9.84801, shape=46.446, scale=.186054, lower.tail=F)
    #    message("Eigenvalue: ", eval, ", median eval: ", qWishartMax(.5, ndf=n, pdim=pdim), 
    #        ", TW stat: ", tws[i], ", median TW: ", qtw(.5), 
    #        ", Gamma pval: ", pvals[i], ", TW pval: ", ptw(tws[i], lower.tail=F))
  }    
  result = list(tws=tws,pvals=pvals)
  return (result)
}

#
# Combine p-values using Fisher's method (unweighted) or weighted Z-method (weighted).
#  
#   p.values: vector or matrix of p-values; one row per variable, one column per p-value
#   method: either "weighted.Z" or "fisher"
#   weights: vector of weights to use with weighted Z-method
#
combinePValues = function(p.values, method="weighted.Z", weights=NA) {    
  if (!(method == "fisher" | method == "weighted.Z")) {
    stop("Method must be either 'fisher' or 'weighted.Z'")
  }  
  if (is.na(weights) & method == "weighted.Z") {
    stop("Weights must be specified for weighted Z-method!")
  }
  
  # Compute the p-values using either Fisher's method or the weighted Z-method
  combined.pvalues = rep(.5, nrow(p.values))
  names(combined.pvalues) = rownames(p.values)  
  for (i in 1:nrow(p.values)) {    
    p = p.values[i,]
    if (method == "fisher") {
      combined.pvalues[i] = pchisq(-2 * sum(log(p)),df=2*length(p),lower.tail=FALSE)
    } else {
      Zi = qnorm(1-p)
      Z  = sum(weights*Zi)/sqrt(sum(weights^2))
      combined.pvalues[i] = 1-pnorm(Z)
    }
  }
  
  return (combined.pvalues)
}

