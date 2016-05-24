######################################################################
#
# expected.deviance.R
#
# copyright (c) 2011-05-30, Katalin Csillery, Olivier Francois and
# Michael GB Blum
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
# Part of the R/abc package
# Contains: expected.deviance
#
###################################################################### 

expected.deviance <- function(target, postsumstat, kernel = "gaussian", subset=NULL, print=TRUE){
    
    if(missing(target)) stop("'target' is missing with no default", call.=F)
    if(missing(postsumstat)) stop("'postsumstat' is missing with no default", call.=F)
    if(!is.matrix(postsumstat) && !is.data.frame(postsumstat) && !is.vector(postsumstat)) stop("'postsumstat' has to be a matrix, data.frame or vector.", call.=F)
    if(!any(kernel == c("gaussian", "epanechnikov", "rectangular", "triangular", "biweight", "cosine"))){
        kernel <- "gaussian"
        warning("Kernel is incorrectly defined. Setting to default (gaussian)", call.=F, immediate.=T)
    }
    if(is.data.frame(postsumstat)) postsumstat <- as.matrix(postsumstat)
    if(is.list(target)) target <- unlist(target)
    if(is.vector(postsumstat)) postsumstat <- matrix(postsumstat, ncol=1)
    if(length(target)!=dim(postsumstat)[2]) stop("Number of summary statistics in 'target' has to be the same as in 'postsumstat'.", call.=F)
    
    ## stop if zero var in postsumstat
    ## #########################
    nss <- length(postsumstat[1,])
    numsim <- length(postsumstat[,1])
    cond1 <- !any(as.logical(apply(postsumstat, 2, function(x) length(unique(x))-1)))
    if(cond1) stop("Zero variance in the summary statistics.", call.=F)

    ## postsumstat values that are to be excluded
    ## #######################################################
    gwt <- rep(TRUE,length(postsumstat[,1]))
    gwt[attributes(na.omit(postsumstat))$na.action] <- FALSE
    if(is.null(subset)) subset <- rep(TRUE,length(postsumstat[,1]))
    gwt <- as.logical(gwt*subset)
    
    ## scale 
    ## #######
    scaled.postsumstat <- postsumstat
    for(j in 1:nss) scaled.postsumstat[,j] <- normalise(postsumstat[,j],postsumstat[,j][gwt])
    for(j in 1:nss) target[j] <- normalise(target[j],postsumstat[,j][gwt])

    ## calculate euclidean distance
    ## ############################
    sum1 <- 0
    for(j in 1:nss){
        sum1 <- sum1 + (scaled.postsumstat[,j]-target[j])^2
    }
    dist <- sqrt(sum1)

    ## epsilon
    epsilon <- max(dist)
    dist <- dist/epsilon

    ## apply kernel
    ## ##############
    if(kernel == "gaussian") dist <- exp(-.5*dist^2)/sqrt(2*pi)
    if(kernel == "epanechnikov") dist <- 1 - (dist)^2
    if(kernel == "rectangular") dist <- dist ## uniform
    if(kernel == "triangular") dist <- 1 - abs(dist)
    if(kernel == "biweight") dist <- (1 - (dist)^2)^2
    if(kernel == "cosine") dist <- cos(pi/2*dist)
    dist <- dist/epsilon
    dist[dist==0] <- NA
    if(print) cat(sum(is.na(dist)), "/", length(dist)," replicated data were excluded.\n", sep="")
    dist[is.na(dist)] <- 4e-44
    dist <- log(dist)
    
    ## calculate expected.deviance
    ## ###################
    dev <- - 2*mean(dist)	
    return(list(expected.deviance=dev, dist=dist))
}
