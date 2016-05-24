### pvt.cppowscor.R  (2012-01-21)
###
###    Efficient computation of crossprod(R^alpha, y)
###
### Copyright 2011-2012 A. Pedro Duarte Silva, Verena Zuber, and Korbinian Strimmer
###
###
### This file is part of the `corpcor' library for R and related languages.
### It is made available under the terms of the GNU General Public
### License, version 3, or at your option, any later version,
### incorporated herein by reference.
### 
### This program is distributed in the hope that it will be
### useful, but WITHOUT ANY WARRANTY; without even the implied
### warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
### PURPOSE.  See the GNU General Public License for more
### details.
### 
### You should have received a copy of the GNU General Public
### License along with this program; if not, write to the Free
### Software Foundation, Inc., 59 Temple Place - Suite 330, Boston,
### MA 02111-1307, USA





##### internal functions ######

# this procedure exploits a special identity to efficiently
# compute the crossprod of matrix power of the correlation shrinkage estimator with y

# computes R_shrink^alpha %*% y
pvt.cppowscor = function(x, y, alpha, lambda, w, verbose)
{
  # determine correlation shrinkage intensity
  if (missing(lambda))
  {
    lambda = estimate.lambda(x, w, verbose)
    lambda.estimated=TRUE
  }
  else
  {
    if (lambda < 0) lambda = 0
    if (lambda > 1) lambda = 1
    if (verbose)
    {
      cat(paste("Specified shrinkage intensity lambda (correlation matrix):", round(lambda, 4), "\n"))     
    }
    lambda.estimated=FALSE
  }

  #####
 
  n = nrow(x)  
  w = pvt.check.w(w, n)
 
  # standardize input matrix by standard deviations
  xs = wt.scale(x, w, center=TRUE, scale=TRUE) # standardize data matrix

  # bias correction factor
  h1 = 1/(1-sum(w*w))   # for w=1/n this equals the usual h1=n/(n-1)

  p = ncol(xs)
    
  if (lambda == 1 | alpha == 0) # result in both cases R is the identity matrix
  {
      cp.powr = y    # return y
  }
  else
  {
    # number of zero-variance variables
    zeros = (attr(xs, "scaled:scale")==0.0)
    
    svdxs = fast.svd(xs)
    m = length(svdxs$d)  # rank of xs
               
    UTWU = t(svdxs$u) %*% sweep(svdxs$u, 1, w, "*") #  t(U) %*% diag(w) %*% U
    C = sweep(sweep(UTWU, 1, svdxs$d, "*"), 2, svdxs$d, "*") # D %*% UTWU %*% D
    C = (1-lambda) * h1 * C

    C = (C + t(C))/2  # symmetrize for numerical reasons (mpower() checks symmetry)
    
    # note: C is of size m x m, and diagonal if w=1/n
         
    if (lambda==0.0) # use eigenvalue decomposition computing the matrix power
    {
      if (m < p-sum(zeros)) 
        warning(paste("Estimated correlation matrix doesn't have full rank",
          "- pseudoinverse used for inversion."), call. = FALSE)    
       
      cp.powr =  svdxs$v %*%  (mpower(C, alpha)  %*% crossprod( svdxs$v, y))
    }
    else # use a special identity for computing the matrix power
    {
      F = diag(m) - mpower(C/lambda + diag(m), alpha)
      cp.powr = (y - svdxs$v %*% (F %*% crossprod(svdxs$v, y) ))*(lambda)^alpha 
    }
 
    # set all diagonal entries in R_shrink corresponding to zero-variance variables to 1
    cp.powr[zeros,] = y[zeros,]
  }
  rownames(cp.powr) = colnames(xs)
  colnames(cp.powr) = colnames(y)  
  rm(xs)

  attr(cp.powr, "lambda") = lambda
  attr(cp.powr, "lambda.estimated") = lambda.estimated
  attr(cp.powr, "class") = "shrinkage"

  if (verbose) cat("\n")
  
  return( cp.powr )
}




