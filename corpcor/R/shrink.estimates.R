### shrink.estimates.R  (2012-01-21)
###
###    Shrinkage Estimation of Variance Vector, Correlation Matrix,
###    and Covariance Matrix 
###
### Copyright 2005-12 Juliane Schaefer, Rainer Opgen-Rhein, 
###        Verena Zuber, A. Pedro Duarte Silva, and Korbinian Strimmer
###
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



# power of the shrinkage correlation matrix
powcor.shrink = function(x, alpha, lambda, w, verbose=TRUE)
{
  if (missing(alpha)) stop("Please specify the exponent alpha!")

   x = as.matrix(x)
   
   # matrix power of shrinkage correlation
   powr = pvt.powscor(x=x, alpha=alpha, lambda=lambda, w=w, verbose=verbose)
   
   return(powr)
}


# correlation
cor.shrink = function(x, lambda, w, verbose=TRUE)
{
   return ( powcor.shrink(x=x, alpha=1, lambda=lambda, w=w, verbose=verbose) )
}


# inverse correlation
invcor.shrink = function(x, lambda, w, verbose=TRUE)
{
   return ( powcor.shrink(x=x, alpha=-1, lambda=lambda, w=w, verbose=verbose) )
}


# variances
var.shrink = function(x, lambda.var, w, verbose=TRUE)
{
  x = as.matrix(x)
  
  # shrinkage variance 
  sv = pvt.svar(x=x, lambda.var=lambda.var, w=w, verbose=verbose)
  
  return(sv)
}


# covariance
cov.shrink = function(x, lambda, lambda.var, w, verbose=TRUE)
{   
   x = as.matrix(x)
   
   # shrinkage scale factors
   sc = sqrt( pvt.svar(x=x, lambda.var=lambda.var, w=w, verbose=verbose) )

   # shrinkage correlation
   c = pvt.powscor(x=x, alpha=1, lambda=lambda, w=w, verbose=verbose)
   
   # shrinkage covariance 
   if (is.null(dim(c)))
     c = c*sc*sc
   else
     c = sweep(sweep(c, 1, sc, "*"), 2, sc, "*")

   attr(c, "lambda.var") = attr(sc, "lambda.var")
   attr(c, "lambda.var.estimated") = attr(sc, "lambda.var.estimated")
                    
   return(c)
}


# precision matrix (inverse covariance)
invcov.shrink = function(x, lambda, lambda.var, w, verbose=TRUE)
{   
   x = as.matrix(x)

   # shrinkage scale factors
   sc = sqrt( pvt.svar(x=x, lambda.var=lambda.var, w=w, verbose=verbose) )
        
   # inverse shrinkage correlation
   invc = pvt.powscor(x=x, alpha=-1, lambda=lambda, w=w, verbose=verbose)
   
   # inverse shrinkage covariance 
   if (is.null(dim(invc)))
     invc = invc/sc/sc
   else
     invc = sweep(sweep(invc, 1, 1/sc, "*"), 2, 1/sc, "*")
   
   attr(invc, "lambda.var") = attr(sc, "lambda.var")
   attr(invc, "lambda.var.estimated") = attr(sc, "lambda.var.estimated")
   
   return(invc)
}


# computes R_shrink^alpha %*% y
crossprod.powcor.shrink = function(x, y, alpha, lambda, w, verbose=TRUE)
{
  if (missing(alpha)) stop("Please specify the exponent alpha!")

   x = as.matrix(x)
   y = as.matrix(y)
   p = ncol(x)
   if (nrow(y) != p) stop("Matrix y must have ", p, " rows!")
   
   # crossprod of matrix power of shrinkage correlation with y
   cp.powr = pvt.cppowscor(x=x, y=y, alpha=alpha, lambda=lambda, w=w, verbose=verbose)
   
   return(cp.powr)
}


