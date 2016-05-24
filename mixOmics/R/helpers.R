# Copyright (C) 2015
# Kim-Anh Le Cao, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD
# Florian Rohart, Australian Institute for Bioengineering and Nanotechnology, University of Queensland, Brisbane, QLD.
# Benoit Gautier, The University of Queensland, The University of Queensland Diamantina Institute, Translational Research Institute, Brisbane, QLD

# This program is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 2
# of the License, or (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, write to the Free Software
# Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA  02111-1307, USA.



# ----------------------------------------------------------------
# Internal helpers functions to run the GCCA module, not to be used externally
# ---------------------------------------------------------------- 
# Some of these functions have been borrowed from the RGCCA package, as indicated below

# --------------------------------------
# soft_thresholding
# --------------------------------------
soft_thresholding_L1=function(x,nx){
  
  #selection on a (loadings.X). modified on 19/02/15 to make sure that a!=0
  if(nx!=0)
  {
    absa = abs(x)
    if(any(rank(absa, ties.method = "max") <= nx)) {
      x = ifelse(rank(absa, ties.method = "max") <= nx, 0,
                 sign(x) * (absa - max(absa[rank(absa, ties.method = "max") <= nx])))
    }
  }
  
  x
}

# ----------------------------------------------------------------------------------------------------------
# soft.threshold() - soft-thresholds a vector such that the L1-norm constraint is satisfied.
# ----------------------------------------------------------------------------------------------------------
soft.threshold = function (x, sumabs = 1)
  return(soft(x, BinarySearch(x, sumabs)))

BinarySearch = function(argu,sumabs){
  if(norm2(argu)==0 || sum(abs(argu/norm2(argu)))<=sumabs) return(0)
  lam1 <- 0
  lam2 <- max(abs(argu))-1e-5
  iter <- 1
  while(iter < 500){
    su <- soft(argu,(lam1+lam2)/2)
    if(sum(abs(su/norm2(su)))<sumabs){
      lam2 <- (lam1+lam2)/2
    } else {
      lam1 <- (lam1+lam2)/2
    }
    if((lam2-lam1)<1e-10) return((lam1+lam2)/2)
    iter <- iter+1
  }
  warning("Didn't quite converge")
  return((lam1+lam2)/2)
}

soft = function(x,d) return(sign(x)*pmax(0, abs(x)-d))

norm2 = function(vec){
  a <- sqrt(sum(vec^2))
  if(a==0) a <- .05
  return(a)
}


# --------------------------------------
# sparsity function
# --------------------------------------
sparsity=function(loadings.A, keepA, penalty){
  
  if(!is.null(keepA)) {
    nx=length(loadings.A) - keepA
    loadings.A = soft_thresholding_L1(loadings.A,nx)
  } else if (!is.null(penalty)){
    loadings.A = soft.threshold(loadings.A, penalty)
  }
  return(loadings.A)
}



# --------------------------------------
# l2.norm
# --------------------------------------
l2.norm=function(x)
{
  if(!is.vector(x)) stop("x has to be a vector")
  out=x/drop(sqrt(crossprod(x)))
}


#############################################################################################################
# Functions acquired from RGCCA R-library
#############################################################################################################
# ----------------------------------------------------------------------------------------------------------
# cov2() - Compute biased and unbiased covariance and variance estimates
# ----------------------------------------------------------------------------------------------------------
cov2 = function (x, y = NULL, bias = TRUE) {
  n = NROW(x)
  if (is.null(y)) {
    x = as.matrix(x)
    if (bias) {
      C = ((n - 1)/n) * cov(x, use = "pairwise.complete.obs")
    } else {
      C = cov(x, use = "pairwise.complete.obs")
    }
  } else {
    if (bias) {
      C = ((n - 1)/n) * cov(x, y, use = "pairwise.complete.obs")
    } else {
      C = cov(x, y, use = "pairwise.complete.obs")
    }
  }
  return(C)
}

# ----------------------------------------------------------------------------------------------------------
# scale2() - scales columns (variables) of each datasets in list A
# ----------------------------------------------------------------------------------------------------------
scale2 = function (A, center = TRUE, scale = TRUE, bias = TRUE) {
  ### Start: add row.names to the output
  rownames = NULL
  if (!is.null(row.names(A)))
    rownames = row.names(A)
  ### End: add row.names to the output
  
  if (center == TRUE & scale == TRUE) {
    A = scale(A, center = TRUE, scale = FALSE)
    std = sqrt(apply(A, 2, function(x) cov2(x, bias = bias)))
    #if (any(std == 0)) {
    #    sprintf("there were %d constant variables", sum(std ==
    #    0))
    #    std[std == 0] = 1
    #}
    A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    ### return(A) Wait end script
  }
  if (center == TRUE & scale == FALSE) {
    A = scale(A, center = TRUE, scale = FALSE)
    ### return(A) Wait end script
  }
  if (center == FALSE & scale == TRUE) {
    std = apply(A, 2, function(x) cov2(x, bias = bias))
    A = A/matrix(rep(std, NROW(A)), NROW(A), NCOL(A), byrow = TRUE)
    attr(A, "scaled:scale") = std
    ### return(A) Wait end script
  }
  ### Start: add row.names to the output
  if (!is.null(rownames))
    row.names(A) = rownames
  return(A)
  ### End: add row.names to the output
}

# ----------------------------------------------------------------------------------------------------------
# initsvd() - performs SVD on matrix X
# ----------------------------------------------------------------------------------------------------------
initsvd = function (X) {
  n = NROW(X)
  p = NCOL(X)
  ifelse(n >= p, return(svd(X, nu = 0, nv = 1)$v), return(svd(X, nu = 1, nv = 0)$u))
}

# ----------------------------------------------------------------------------------------------------------
# miscrossprod() - Compute cross-product between vectors x and y
# ----------------------------------------------------------------------------------------------------------
miscrossprod = function (x, y) {
  d.p = sum(drop(x) * drop(y), na.rm = TRUE)
  #d.p = as.vector(d.p)/norm2(d.p)     ## change made
  return(d.p)
}


# ----------------------------------------------------------------------------------------------------------
# defl.select() - computes residual matrices
# ----------------------------------------------------------------------------------------------------------
defl.select = function(yy, rr, nncomp, nn, nbloc, indY = NULL, mode = "canonical", aa = NULL) { ### Start: Add new parameter for estimation classic mode
  resdefl <- NULL
  pdefl <- NULL
  for (q in 1 : nbloc) {
    ### Start: insertion of new deflations (See La regression PLS Theorie et pratique (page 139))
    if ( nn <= nncomp[q] ) {
      if ((mode == "canonical") || (q != indY)) {
        defltmp <- deflation(rr[[q]], yy[ , q])
        resdefl[[q]] <- defltmp$R
        pdefl[[q]]   <- defltmp$p
      } else if (mode == "classic") {
        resdefl[[q]] <- Reduce("+", lapply(c(1:nbloc)[-q], function(x) {rr[[q]] - yy[ ,x]%*%t(aa[[q]])}))/(nbloc-1)
        pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
      } else if (mode == "invariant") {
        resdefl[[q]] <- rr[[q]]
        pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
      } else if (mode == "regression") {
        resdefl[[q]] <- Reduce("+", lapply(c(1:nbloc)[-q], function(x) {deflation(rr[[q]],yy[, x])$R}))/(nbloc-1)
        pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
      }
      ### End: insertion of new deflations (See La regression PLS Theorie et pratique (page 139))
    } else {
      resdefl[[q]] <- rr[[q]]
      pdefl[[q]]   <-  rep(0,NCOL(rr[[q]]))
    }
  }
  return(list(resdefl=resdefl,pdefl=pdefl))
}


deflation <- function(X, y){
  # Computation of the residual matrix R
  # Computation of the vector p.
  if (any(is.na(X)) || any(is.na(y))) {
    p <- apply(t(X),1,miscrossprod,y)/as.vector(crossprod(y))
  } else {
    p <- t(X)%*%y/as.vector(crossprod(y))
  }
  
  R <- X - y%*%t(p)
  return(list(p=p,R=R))
}
