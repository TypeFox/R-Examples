### carscore.R  (2012-01-21)
###
###    Estimate CAR scores and marginal correlations
###
### Copyright 2010-2012 Verena Zuber and Korbinian Strimmer
###
###
### This file is part of the `care' library for R and related languages.
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




# estimate car scores 
# (and marginal correlations if diagonal=TRUE)

carscore = function(Xtrain, Ytrain, lambda, diagonal=FALSE, verbose=TRUE)
{
  n = dim(Xtrain)[1]
  p = dim(Xtrain)[2]

  #####################################

  if( missing(lambda) )
  {
    lambda.estimated = TRUE
    # regularize the joint correlation matrix  Ytrain and Xtrain combined
    lambda = estimate.lambda( cbind(Ytrain,Xtrain), verbose=verbose )
  }
  else
  {
    lambda.estimated = FALSE
  }
  omega = (1-lambda)*cor(Xtrain, Ytrain) # marginal correlations
  
  if (diagonal==FALSE)
  {
      # car score
      omega = crossprod.powcor.shrink(Xtrain, omega, alpha=-1/2, lambda=lambda, verbose=FALSE)
  }

  omega = as.vector(omega)
  names(omega) = colnames(Xtrain)
  #if(lambda > 0)
  #{
    class(omega) = "shrinkage"
    attr(omega, "lambda") = lambda
    attr(omega, "lambda.estimated") = lambda.estimated
  #}

  return( omega )
}

