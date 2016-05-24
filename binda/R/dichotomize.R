### dichotomize.R  (2015-02-28)
###
###    Dichotomize Continuous Labeled Data 
###
### Copyright 2013-15  Sebastian Gibb and Korbinian Strimmer
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


# dichotomize using given threshold

dichotomize = function(X, thresh)
{
  if(length(thresh)==1) 
    thresh = rep(thresh, ncol(X))
  if(length(thresh) != ncol(X)) 
    stop("Number of specified thresholds not identical with the number of variables (columns)!")


  ## compare features (columns) with thresholds
  ## (we need rowwise comparison but R supports only columnwise comparison
  ##  => double t() )
  m = t(t(X) >= thresh) * 1L

  ## set attributes
  attr(m, "thresh") = thresh

  return(m)
}


# optimize threshold

optimizeThreshold = function(X, L, lambda.freqs, verbose=FALSE)
{
  L = factor(L) # make sure L is a factor
  d = ncol(X) # dimension, number of variables, number of columns

  if (verbose) reportDetails(X, L)

  ## calculate class frequencies
  freqs = getClassFreqs(L, lambda.freqs=lambda.freqs, verbose=verbose)

  ## calculate grid of possible thresholds
  ##  rows: thresholds 1:length(breaks); columns: features
 
  # we simply use all possible thresholds
  grid = apply(X, 2, sort ) # sort to get the smallest cutoff value
  grid = rbind(grid, grid[nrow(X),]+1)
  breaks = nrow(grid)

  ## create score matrix
  ##  rows: thresholds 1:length(breaks); columns: variables
  scr = matrix(0L, nrow=breaks, ncol=d)

  ## loop through thresholds (columns)
  for (i in 1:breaks)
  {
    ## create binary matrix
    ## compare features (columns) with thresholds
    ## (we need rowwise comparison but R supports only columnwise comparison
    ##  => double t() )
    bm = t(t(X) >= grid[i, ]) * 1L

    ## mu matrix
    mu = getClassMeans(bm, L)
 
    scr[i, ] = rankingScore(mu, freqs)
  }

  ## find thesholds with maximal scores
  idx = cbind(max.col(t(scr), ties.method="first"), 1:d)
  thr = grid[idx]

  names(thr) = colnames(X)

  return(thr)
}



