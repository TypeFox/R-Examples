
# This library is free software; you can redistribute it and/or
# modify it under the terms of the GNU Library General Public
# License as published by the Free Software Foundation; either
# version 2 of the License, or (at your option) any later version.
#
# This library is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Library General Public License for more details.
#
# You should have received a copy of the GNU Library General
# Public License along with this library; if not, write to the
# Free Foundation, Inc., 59 Temple Place, Suite 330, Boston,
# MA  02111-1307  USA


################################################################################
# FUNCTION:             DESCRIPTION:
#  .cov.shrink           Builtin from Package 'corpcor'
#  .cor.shrink
#  .varcov
#  .cov.bagged           Builtin from Package 'corpcor'
#  .cor.bagged
#  .bag.fun
#  .robust.cov.boot
#  .sm2vec
#  .smindexes
#  .vec2sm
################################################################################


# Rmetrics:
#   Note that corpcor is not available on Debian as of 2009-04-28. 
#   To run these functions under Debian/Rmetrics we have them    
#   implemented here as a builtin.
#   We also made modifications for tailored usage with Rmetrics. 


# Package: corpcor
# Version: 1.1.2
# Date: 2005-12-12
# Title: Efficient Estimation of Covariance and (Partial) Correlation
# Author: Juliane Schaefer <schaefer@stat.uni-muenchen.de> and
#   Korbinian Strimmer <korbinian.strimmer@lmu.de>.
# Maintainer: Korbinian Strimmer <korbinian.strimmer@lmu.de>
# Depends: R (>= 2.0.0)
# Suggests: 
# Description: This package implements a shrinkage estimator to allow
#   the efficient inference of large-scale covariance matrices 
#   from small sample data.  The resulting estimates are always
#   positive definite, more accurate than the empirical estimate,
#   well conditioned, computationally inexpensive, and require
#   only little a priori modeling.  The package also contains
#   similar functions for inferring correlations and partial
#   correlations.  In addition, it provides functions for fast svd 
#   computation, for computing the pseuoinverse, and 
#   for checking the rank and positive definiteness of a matrix.
# License: GPL version 2 or newer
# URL: http://www.statistik.lmu.de/~strimmer/software/corpcor/
# Packaged: Mon Dec 12 13:07:22 2005; strimmer


# ------------------------------------------------------------------------------


.cov.shrink <- 
function(x, lambda, verbose = FALSE)
{
   x = as.matrix(x)

   # Shrinkage correlation coefficients
   R.star <- .cor.shrink(x, lambda = lambda, verbose=verbose)

   # Unbiased empirical variances
   V = apply(x, 2, var)
   
   resid.sd = sqrt(V)
   ans <- sweep(sweep(R.star, 1, resid.sd, "*"), 2, resid.sd, "*") 
     
   # Return Value:
   ans
}


# ------------------------------------------------------------------------------


.cor.shrink <- 
function(x, lambda, verbose = FALSE)
{
    # Standardize data (and turn x into a matrix)
    sx <- scale(x)  
    
    p = dim(sx)[2]
    if(p == 1) return( as.matrix(1) ) 
    
    # Estimate variance of empirical correlation coefficients 
    vc = .varcov(sx, type = "unbiased", verbose)
    
    # Find optimal lambda:
    if(missing(lambda)) {   
        offdiagsum.rij.2 = sum(vc$S[lower.tri(vc$S)]^2)
        offdiagsum.v.rij = sum(vc$var.S[lower.tri(vc$var.S)])     
        lambda = offdiagsum.v.rij/offdiagsum.rij.2
        if(verbose) cat(paste("Estimated shrinkage intensity lambda: ",
            round(lambda,4), "\n"))
    }
    if(lambda > 1) {
        warning(paste("Overshrinking: lambda set to 1 (allowed range: 0-1)"))
        lambda = 1  
    } else if(lambda < 0) {
        warning(paste("Undershrinking: lambda set to 0 (allowed range: 0-1)"))
        lambda = 0  
    }
 
    # construct shrinkage estimator
    R.star = (1-lambda) * vc$S
    diag(R.star) = rep(1, p)
    attr(R.star, "lambda") = lambda
  
    # Return Value:
    R.star
}


# ------------------------------------------------------------------------------


.varcov <-  
function(x, type = c("unbiased", "ML"), verbose = FALSE)
{
    # Details:
    #   compute the empirical covariance matrix S=cov(x) given a data 
    #   matrix x as well as the *variances* associated with the individual 
    #   entries S[i,j]

    x = as.matrix(x)     
    n = dim(x)[1]
    p = dim(x)[2]
         
    # Weights for the "unbiased" and "ML" cases
    type = match.arg(type)
    if(type == "unbiased") {
        h1 = 1/(n-1)
        h2 = n/(n-1)/(n-1) }    
    if(type == "ML") {
      h1 = 1/n
      h2 = (n-1)/n/n
    }
 
    s = matrix(NA, ncol = p, nrow = p)   
    vs = matrix(NA, ncol = p, nrow = p)
    xc = scale(x, scale=FALSE) # center the data
    
    # Diagonal elements:
    for (i in 1:p) {
        zii = xc[,i]^2
        s[i, i] = sum(zii)*h1
        vs[i, i] = var(zii)*h2
    }
    if(p == 1) 
        return(list(S = s, var.S = vs))
    if(verbose && p > 50)
        cat(paste("Computing ... wait for", p, "dots (50 per row):\n")) 
    
    # Off-diagonal elements
    for (i in 1:(p-1)) {
        if(verbose && p > 50) {
            cat(".")
            if(i %% 50 == 0) cat(paste(" ", i, "\n"))
        }
        for (j in (i+1):p) {
            zij = xc[,i]*xc[, j] 
            s[i, j] = sum(zij)*h1
            s[j, i] = s[i,j]
            vs[i, j] = var(zij)*h2
            vs[j, i] = vs[i, j]   
        }
    }
    if(verbose && p > 50) cat(paste(". ", i+1, "\n"))

    # Return Value:
    return(list(S = s, var.S = vs))
}


################################################################################
# cov.bagged.R  (2004-03-15)
#   Variance reduced estimators of cov, cor, and pcor
#       using bootstrap aggregation ("bagging")
#   Copyright 2003-04 Juliane Schaefer and Korbinian Strimmer
# Package: corpcor
# Version: 1.1.2
# Date: 2005-12-12
# Title: Efficient Estimation of Covariance and (Partial) Correlation
# Author: Juliane Schaefer <schaefer@stat.uni-muenchen.de> and
#         Korbinian Strimmer <korbinian.strimmer@lmu.de>.
# Maintainer: Korbinian Strimmer <korbinian.strimmer@lmu.de>
# Depends: R (>= 2.0.0)
# Suggests: 
# Description: This package implements a shrinkage estimator to allow
#   the efficient inference of large-scale covariance matrices 
#   from small sample data.  The resulting estimates are always
#   positive definite, more accurate than the empirical estimate,
#   well conditioned, computationally inexpensive, and require
#   only little a priori modeling.  The package also contains
#   similar functions for inferring correlations and partial
#   correlations.  In addition, it provides functions for fast svd 
#   computation, for computing the pseuoinverse, and 
#   for checking the rank and positive definiteness of a matrix.
# License: GPL version 2 or newer
# URL: http://www.statistik.lmu.de/~strimmer/software/corpcor/
# Packaged: Mon Dec 12 13:07:22 2005; strimmer


.cov.bagged <-  
function(x, R = 1000, ...)
{
    vec.out = .bag.fun(cov, x, R = R, diag = TRUE, ...)
    mat.out = .vec2sm(vec.out, diag = TRUE)
  
    # Return Value:
    mat.out
}


# ------------------------------------------------------------------------------


.cor.bagged <-  
function(x, R = 1000, ...)
{
    vec.out = .bag.fun(cor, x, R = R, diag = FALSE, ...)
    mat.out = .vec2sm(vec.out, diag = FALSE)
    
    # Fill diagonal with 1
    diag(mat.out) = rep(1, dim(mat.out)[1]) 
  
    # Return Value:
    mat.out
}


# ------------------------------------------------------------------------------


.bag.fun <-  
function(fun, data, R, diag, ...)
{
    # Number of variables 
    p = dim(data)[2]
  
    # Index vector for lower triangle
    lo = lower.tri(matrix(NA, nrow=p, ncol=p), diag=diag)

    # bootstrap function
    .bootFun = function(data, i) {
        vec = as.vector( fun(data[i,], ...)[lo] )
        # if we get NAs flag result as being erroneous
        if(sum(is.na(vec)) > 0) class(vec) = "try-error"
        return( vec )
    }   
     
    # Bag variable 
    boot.out = .robust.cov.boot(data = data, statistic = .bootFun, R = R)
    bag = apply( boot.out$t, 2, mean)
    
    # Return Value:
    bag
}


# ------------------------------------------------------------------------------


.robust.cov.boot <-  
function(data, statistic, R)
{
    # Description:
    #   Simple bootstrap function (robust against errors)
    
    idx = 1:dim(data)[1]
  
    # Determine dimension of statistic
    repeat {
        bx = sample(idx, replace = TRUE)
        val = try(statistic(data, bx))  
        if(class(val) != "try-error") break
    }
    dim.statistic = length(val)
    output = matrix(nrow = R, ncol = dim.statistic)
  
    replicate.count = 0
    error.count = 0
    
    while (replicate.count < R) {
        bx = sample(idx, replace=TRUE)
        val = try(statistic(data, bx)) 
        # if we get a numerical error we simply repeat the draw ..
        if(class(val) == "try-error") {
            error.count = error.count+1   
            if(error.count > R) 
                stop("Too many errors encountered during the bootstrap.")
        } else {
            replicate.count = replicate.count+1
            output[replicate.count, ] = val
        }
    }
    
    if(error.count > 0)  {
        warning(paste(error.count, "out of", R,
            "bootstrap samples were repeated due to errors."))
    }
    
    # Result:
    ans = list(t = output)
            
    # Return Value:
    ans
} 


################################################################################
# smtools.R  (2004-01-15)
#   Convert symmetric matrix to vector and back
#   Copyright 2003-04 Korbinian Strimmer
#
# This file is part of the 'corpcor' library for R and related languages.
# It is made available under the terms of the GNU General Public
# License, version 2, or at your option, any later version,
# incorporated herein by reference.
# 
# This program is distributed in the hope that it will be
# useful, but WITHOUT ANY WARRANTY; without even the implied
# warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
# PURPOSE.  See the GNU General Public License for more
# details.
# 
# You should have received a copy of the GNU General Public
# License along with this program; if not, write to the Free
# Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
# MA 02111-1307, USA


.sm2vec = 
function(m, diag = FALSE)
{
    # Description:
    #   Convert symmetric matrix to vector
    
    ans = as.vector(m[lower.tri(m, diag)])
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.smindexes <-  
function(m, diag = FALSE)
{
    # Descriiption:
    #   Corresponding indexes
      
    m.dim = length(diag(m))
 
    if(diag == TRUE) {
        num.entries = m.dim*(m.dim+1)/2
    } else {
        num.entries = m.dim*(m.dim-1)/2
    }   
    
    index1 = rep(NA, num.entries )
    index2 = rep(NA, num.entries )

    if(diag == TRUE) {
        delta = 0
    } else {
        delta = 1
    }

    z = 1
    for (i in 1:(m.dim-delta)) {
        for (j in (i+delta):m.dim) {
            index1[z] = i
            index2[z] = j
            z = z+1
        }
    }
    
    ans = cbind(index1, index2) 
    
    # Return Value:
    ans
}


# ------------------------------------------------------------------------------


.vec2sm <-  
function(vec, diag = FALSE, order = NULL)
{
    # Description:
    #   Convert vector to symmetric matrix
    
    # Note:
    #   If diag=FALSE then the diagonal will consist of NAs
    
    # dimension of matrix
    n = (sqrt(1+8*length(vec))+1)/2
    if(diag == TRUE) n = n-1
    if( ceiling(n) != floor(n) )
        stop("Length of vector incompatible with symmetric matrix")
       
    # fill lower triangle of matrix     
    m = matrix(NA, nrow = n, ncol = n)
    lo = lower.tri(m, diag)
    if(is.null(order)) {
        m[lo] = vec
    } else {
        # sort vector according to order
        vec.in.order = rep(NA, length(order))
        vec.in.order[order] = vec
        m[lo] = vec.in.order
    }
  
    # symmetrize
    for (i in 1:(n-1)) {
        for (j in (i+1):n) {
            m[i, j] = m[j, i]
        }
    }   
  
    # Return Value:
    m
}


################################################################################

