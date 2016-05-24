### chances.R  (2014-04-08)
###
###    Estimate Bernoulli Parameters from Binary Matrix with Class Labels
###
### Copyright 2014  Sebastian Gibb and Korbinian Strimmer
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



chances = function(X, L, lambda.freqs, verbose=TRUE)
{
  if (!is.binaryMatrix(X)) stop("Input matrix must contain only 0s or 1s!")
  if (missing(L)) stop("Class labels are missing!")
  L = factor(L) # make sure L is a factor

  if (verbose) reportDetails(X, L)

  # class samples
  samples = getClassSamples(L)

  # class frequencies
  regularization = rep(NA, 1)
  names(regularization) = c("lambda.freqs")
  freqs=getClassFreqs(L, lambda.freqs=lambda.freqs, verbose=verbose)
  regularization[1] = attr(freqs, "lambda.freqs")
  attr(freqs, "lambda.freqs")=NULL

  # means
  mu = getClassMeans(X, L)

  # pooled mean 
  mu.pooled = mu %*% freqs
  colnames(mu.pooled) = c("(pooled)")

  means = cbind(mu, mu.pooled)

  out = list(samples=samples, regularization=regularization, freqs=freqs, means=means)

  return(out)
}


