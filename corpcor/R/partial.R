### partial.R  (2008-11-14)
###
###    Partial Correlation and Partial Variance
###    
###
### Copyright 2003-2008 Juliane Schaefer and Korbinian Strimmer
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


#
# partial correlation matrix
#
# input: covariance matrix or correlation matrix
# ouput: partial correlation matrix
#
cor2pcor = function(m, tol)
{   
  # invert, then negate off-diagonal entries
  m = -pseudoinverse(m, tol=tol)
  diag(m) = -diag(m)

  # standardize and return  
  return(cov2cor(m))
}


#
# backtransformation to correlation matrix
#
# input: partial correlation matrix
# ouput: correlation matrix
pcor2cor = function(m, tol)
{
  # negate off-diagonal entries, then invert
  m = -m
  diag(m) = -diag(m)
  m = pseudoinverse(m, tol=tol)
  
  # standardize and return 
  return(cov2cor(m))
}



########################################################


# partial correlation
pcor.shrink = function(x, lambda, w, verbose=TRUE)
{
  pc = -invcor.shrink(x, lambda, w, verbose=verbose)
  diag(pc) = -diag(pc)

  spv = 1/diag(pc)  # standardized partial variances (i.e. pvar/var)
  pc = cov2cor(pc)  # partial correlations

  attr(pc, "spv") = spv
 
  return(pc)
}


# partial variances
pvar.shrink = function(x, lambda, lambda.var, w, verbose=TRUE)
{
  prec = invcov.shrink(x, lambda, lambda.var, w, verbose=verbose)
  pvar = 1/diag(prec)

  attr(pvar, "lambda") = attr(prec, "lambda")
  attr(pvar, "lambda.estimated") = attr(prec, "lambda.estimated")
  attr(pvar, "lambda.var") = attr(prec, "lambda.var")
  attr(pvar, "lambda.var.estimated") = attr(prec, "lambda.var.estimated")
    
  return( pvar )
}
