### predict.slm.R  (2014-11-22)
###
###    Prediction from linear model
###
### Copyright 2011-14 Korbinian Strimmer
###
###
### This file is part of the `sda' library for R and related languages.
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




predict.slm = function(object, Xtest, verbose=TRUE, ...)
{
  if ( missing(object) ) {
    stop("An slm object must be supplied.")
  }

  if ( missing(Xtest) ) {
    stop("A test data set must be supplied.")
  }
  
  if (!is.matrix(Xtest)) stop("Test data must be given as matrix!")
  ntest = nrow(Xtest) # number of test samples
  nvtest = ncol(Xtest) # number of of variables in test data set
  ncoeff =  ncol(object$coefficients)-1 # number of coefficients

  if (ncoeff != nvtest)
    stop("Incompatible number of variables in test data set (", nvtest,
          ") and number of coefficients in slm object (", ncoeff, ")", sep="")

  m = length(object$numpred)
  yhat = matrix(0, nrow=ntest, ncol=m)
  colnames(yhat) = names(object$numpred)
  rownames(yhat) = rownames(Xtest)

  for (i in 1:m)
  {
    if (verbose) cat("Prediction uses", object$numpred[i], "variables.\n")

    b = matrix(object$coefficients[i, -1]) 
    b0 = object$coefficients[i, 1]

    yhat[,i] =  b0 + Xtest %*% b 

  }
  return( yhat )
}

