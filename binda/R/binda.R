### binda.R  (2015-02-28)
###
###    Multi-class discriminant analysis with binary predictors
###
### Copyright 2013-2015  Sebastian Gibb and Korbinian Strimmer
###
###
### This file is part of the `binda' library for R and related languages.
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



binda = function(Xtrain, L, lambda.freqs, verbose=TRUE)
{
  if (!is.binaryMatrix(Xtrain)) stop("Training data must be a matrix containing only 0 or 1!")
  if (missing(L)) stop("Class labels are missing!")
  L = factor(L) # make sure L is a factor

  if (verbose) reportDetails(Xtrain, L)

  # class frequencies
  freqs = getClassFreqs(L, lambda.freqs=lambda.freqs, verbose=verbose)
  # store lambda.freqs attr in regularization
  regularization = setNames(attr(freqs, "lambda.freqs"), "lambda.freqs")
  # remove lambda.freqs attribute
  attr(freqs, "lambda.freqs") = NULL

  logfreqs = log(freqs)

  # means
  mu = avoidBoundaries(getClassMeans(Xtrain, L))
  logp1 = log(mu)
  logp0 = log(1-mu)

  out = list(regularization=regularization, logfreqs=logfreqs, logp0=logp0, logp1=logp1)
  class(out) = "binda"

  return( out )
}


predict.binda = function(object, Xtest, verbose = TRUE, ...)
{
  if (missing(object)) {
    stop("A binda fit object must be supplied.")
  }
  if (missing(Xtest)) {
    stop("A new data to predict must be supplied.")
  }
  if (!is.binaryMatrix(Xtest)) stop("Test data must be a matrix containing only 0 or 1!")

  logp0 = object$logp0
  logp1 = object$logp1
  logfreqs = object$logfreqs

  if (ncol(Xtest) != nrow(logp1))
      stop("Different number of predictors in binda object (",
           nrow(logp1), ") and in test data (", ncol(Xtest),
           ")", sep = "")

  if (verbose) cat("Prediction uses", nrow(logp1), "features.\n")

  # compute discriminant function  (rows: test samples, columns: classes)
  probs = matrix( logfreqs + colSums(logp0),
            nrow = nrow(Xtest), ncol = length(logfreqs), byrow=TRUE ) +
            crossprod( t(Xtest), logp1-logp0 )

  # class prediction (with random ties)
  yhat = max.col(probs) # yhat = apply(probs, 1, which.max)

  # convert to probabilities
  probs = exp(probs - max.col.value(probs))  #probs = exp(probs - apply(probs, 1, max))
  probs = zapsmall( probs / rowSums(probs) )

  attr(yhat, "levels") = names(logfreqs)
  class(yhat) = "factor"
  colnames(probs) = names(logfreqs)
  rownames(probs) = rownames(Xtest)

  return( list(class = yhat, posterior = probs) )
}


## private function
max.col.value = function(x)
{
  return( x[cbind(1L:nrow(x), max.col(x, ties.method="first"))] )
}


# ensure that no mean is exactly 0 or 1
avoidBoundaries = function(mu, n=1e9)
{
  lambda = 1/(n + 1)
  target = 1/2

  return( mu * (1-lambda) + lambda*target )
}

