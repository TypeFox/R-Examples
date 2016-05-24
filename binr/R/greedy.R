#-------------------------------------------------------------------------------
#
# Package binr 
#
# Implementation of greedy algorithm. 
# 
# Sergei Izrailev, 2011-2014
#-------------------------------------------------------------------------------
# Copyright 2011-2014 Collective, Inc.
# 
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# 
# http://www.apache.org/licenses/LICENSE-2.0
# 
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
#-------------------------------------------------------------------------------#-------------------------------------------------------------------------------

#' \code{bins.greedy} - Wrapper around \code{bins.greedy.impl}. Goes over the
#' sorted values of \code{x} left to right and fills the bins with the values until
#' they are about the right size.
#' @param x Vector of numbers. 
#' @param nbins Target number of bins.
#' @param minpts Minimum number of points in a bin. Only used if \code{naive = FALSE}.  
#' @param naive When \code{TRUE}, simply calls \code{bins.greedy.impl} with data 
#'        derived from \code{x}. Otherwise, makes an extra step of marking the values that
#'        by themselves take a whole bin to force the algorithm to place these values
#'        in a bin separately.
#' @name bins.greedy
#' @title Greedy binning algorithm.
#' @seealso \code{\link{binr}}, \code{\link{bins}}, \code{\link{bins.quantiles}} \code{\link{bins.optimize}}
#' @export
#' @rdname bins.greedy
bins.greedy <- function(x, nbins, minpts = floor(0.5 * length(x) / nbins), thresh = 0.8, naive = FALSE)
{
   xtbl <- table(x)
   xval <- sort(unique(x))
   nvals <- length(xval)
   xstp <- vector(nvals, mode="logical")
   binsz <- floor(length(x) / nbins)
   if (naive) return( bins.greedy.impl(xval, xtbl, xstp, binsz, nbins, thresh) )
   
   ixdone <- which(xtbl >= binsz)
   ixkeep <- which(xtbl < binsz)
   nbins2 <- nbins - length(ixdone)
   nx2 <- sum(xtbl) - sum(xtbl[ixdone]) # how many points are left
   nbins.new <- min(nbins2, floor(nx2 / minpts) + 1) # make fewer bins if there's not enough points
   binsz.new <- floor(nx2 / nbins.new)
   istp <- ixdone + 1
   istp <- istp[istp < nvals]
   xstp.new <- xstp
   xstp.new[istp] <- TRUE
   #cuts.new <- bins.greedy.impl(xval[ixkeep], xtbl[ixkeep], xstp.new[ixkeep], binsz.new, nbins.new, thresh)
   cuts.new <- bins.greedy.impl(xval, xtbl, xstp, binsz.new, nbins.new, thresh)
   return(cuts.new)
}

#-------------------------------------------------------------------------------

#' \code{bins.greedy.impl} - Implementation of a single-pass binning algorithm that examines sorted data left to right 
#' and builds bins of the target size. The \code{bins.greedy} wrapper around this function provides a less involved interface.  
#' This is not symmetric wrt direction: symmetric distributions may not have symmetric bins if there are multiple points
#' with the same values. If a single value accounts for more than thresh * binsz points, it will be placed in 
#' a new bin.
#' @param xval Sorted unique values of the data set x. This should be the numeric version of \code{names(xtbl)}.
#' @param xtbl Result of a call to \code{table(x)}.
#' @param xstp Stopping points; if \code{xstp[i] == TRUE}, the \code{i}-th value can't be merged to the \code{(i-1)}-th one. 
#'        \code{xstp[1]} value is ignored.
#' @param binsz Target bin size, i.e., the number of points falling into each bin; for example, \code{floor(length(x) / nbins)}
#' @param thresh Threshold fraction of bin size for the greedy algorithm. 
#'        Suppose there's \code{n < binsz} points in the current bin already. 
#'        Also suppose that the next value V is represented by \code{m} points, and \code{m + n > binsz}. 
#'        Then the algorithm will check if \code{m > thresh * binsz}, and if so, will place the value V into a new bin.
#'        If \code{m} is below the threshold, the points having value V are added to the current bin.
#' @param verbose When \code{TRUE}, prints the number of points falling into the bins.
#' @return A list with the following items:
#' \itemize{
#'    \item{binlo}{ - The "low" value falling into the bin.}
#'    \item{binhi}{ - The "high" value falling into the bin.}
#'    \item{binct}{ - The number of points falling into the bin.}
#'    \item{xtbl}{ - The result of a call to \code{table(x)}.}
#'    \item{xval}{ - The sorted unique values of the data points x. Essentially, a numeric version of \code{names(xtbl)}.}
#' }
#' @export
#' @rdname bins.greedy
bins.greedy.impl <- function(xval, xtbl, xstp, binsz, nbins, thresh, verbose=F)
{   
   nvals <- length(xval)
   if (nvals != length(xtbl) || nvals != length(xstp)) stop("bins.greedy: xval, xtbl and xstp lengths don't match")
   xstp[1] = FALSE
   binlo <- vector(nbins, mode="integer")
   binhi <- vector(nbins, mode="integer")
   binct <- vector(nbins, mode="integer")
   startbin <- TRUE
   k <- 1   # current bin index
   s <- 0   # running number of points in the current bin
   for (i in 1:nvals)
   {
      if (xstp[i]) # can't happen on i == 1
      {
         # wrap up the previous bin
         binhi[k] <- xval[i - 1]
         binct[k] <- s
         s <- 0
         k <- k + 1
         startbin <- TRUE         
      }
      s <- s + xtbl[i]
      if (startbin) 
      {
         binlo[k] <- xval[i]
         startbin <- FALSE
      }
      if (s >= binsz)
      {
         # we exceeded the bin size. If there were points in this bin before, check if
         # the added number of points exceeds the tolerance threshold.
         if (s > xtbl[i] && (xtbl[i] >= thresh * binsz))# || (s - xtbl[i] > thresh * binsz))) 
         {
            # not the first item; complete the previous bin and make this a starting bin
            binhi[k] <- xval[i - 1]
            binct[k] <- s - xtbl[i]
            k <- k + 1
            s <- xtbl[i]
            binlo[k] <- xval[i]
            startbin <- FALSE         
         }
         else
         {
            # done for now
            binhi[k] <- xval[i]
            binct[k] <- s
            k <- k + 1
            s <- 0
            startbin <- TRUE         
         }
      } 
   }
   if (!startbin)
   {
      # unfinished bin
      binhi[k] <- xval[nvals]
      binct[k] <- s
   }
   #if (k > nbins) k <- nbins  # we are NOT guaranteed that everything fits into at most nbins bins
   # Here's an example when things go bad: 
   #    x <- rep(1:10, 10); bins.greedy(x, nbins=9, thresh=0.8)$binct
   # In this case, every value exceeds the threshold and is placed in a separate bin, hence we get
   # more bins than what we asked for.
  
   # Reduce k until a valid bin.
   while (is.na(binct[k])) k <- k - 1
   binlo <- binlo[1:k]
   binhi <- binhi[1:k]
   binct <- binct[1:k]
   
   names(binct) <- paste("[", binlo, ", ", binhi, "]", sep="")
   if (verbose) print(binct)
   return(list(binlo=binlo, binhi=binhi, binct=binct, xtbl=xtbl, xval=xval))
}

#-------------------------------------------------------------------------------
