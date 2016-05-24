#-------------------------------------------------------------------------------
#
# Package binr 
#
# Implementation of algorithms minimizing the error of the split. 
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

#' Algorithms minimizing the binning error function \code{bins.merr}.
#' 
#' \code{bins.move} - Compute the best move of a value from one bin to its neighbor
#' @param xval Sorted unique values of the data set x. This should be the numeric version of \code{names(xtbl)}.
#' @param xtbl Result of a call to \code{table(x)}.
#' @param binlo The "low" value falling into the bin.
#' @param binhi The "high" value falling into the bin.
#' @param binct The number of points falling into the bin.
#' @param target.bins Number of bins desired; this is also the max number of bins.
#' @return A list containing the following items (not all of them may be present):
#' \itemize{
#'    \item{binlo}{ - The "low" value falling into the bin.}
#'    \item{binhi}{ - The "high" value falling into the bin.}
#'    \item{binct}{ - The number of points falling into the bin.}
#'    \item{xtbl}{ - The result of a call to \code{table(x)}.}
#'    \item{xval}{ - The sorted unique values of the data points x. Essentially, a numeric version of \code{names(xtbl)}.}
#'    \item{changed}{ - Flag indicating whether the bins have been modified by the function.}
#'    \item{err}{ - Mean square root error between the resulting counts and ideal bins.}
#'    \item{imax}{ - For the move, merge and split operations, the index of the bin with the maximum gain.}
#'    \item{iside}{ - For the move operation, the side of the move: 0 = left, 1 = right.}
#'    \item{gain}{ - Error gain obtained as the result of the function call.}
#' }
#' @name bins.optimize
#' @aliases bins.move bins.split bins.merge bins.move.iter bins.split.iter bins.merge.iter
#' @export
#' @seealso \code{\link{bins}}, \code{\link{binr}}, \code{\link{bins.greedy}}, \code{\link{bins.quantiles}} 
#' @rdname bins.optimize
bins.move <- function(xval, xtbl, binlo, binhi, binct, target.bins, verbose=F)
{
   nbins <- length(binlo)
   if (nbins != length(binhi) || nbins != length(binct)) stop("bins.move: lengths of binlo, binhi, binct do not match")
   chggain <- vector(2*nbins, mode="double")
   chghilo <- vector(nbins, mode="integer")
   bestgain <- -Inf
   bestside <- NA
   bestindx <- -1
   for (i in 1:nbins)
   {
      # take points from this bin if possible and place in the neighboring bin
      # binlo[i] = index of the smallest item in the bin
      # binhi[i] = index of the largest  item in the bin
      nitems <- binhi[i] - binlo[i] + 1
      if (nitems > 1)
      {
         cntlo <- xtbl[binlo[i]]
         cnthi <- xtbl[binhi[i]]
         # try left
         if (i > 1 && binct[i] -  cntlo > binct[i - 1])
         {
            # compute gain when moving left
            gain <- 2 * cntlo * (binct[i] - cntlo - binct[i - 1]) / nbins
            if (gain > bestgain && gain > 0)
            {
               bestgain <- gain
               bestside <- 0
               bestindx <- i
            }
         }
         # try right
         if (i < nbins && binct[i] - cnthi > binct[i + 1])
         {
            # compute gain when moving right
            gain <- 2 * cntlo * (binct[i] - cnthi - binct[i + 1]) / nbins
            if (gain > bestgain && gain > 0)
            {
               bestgain <- gain
               bestside <- 1
               bestindx <- i
            }
         }
      }
   }
   
   changed = FALSE
   if (!is.na(bestside))
   {
      i <- bestindx
      cntlo <- xtbl[binlo[i]]
      cnthi <- xtbl[binhi[i]]
      if (bestside == 1)
      {
         # move right
         binct[i + 1] <- binct[i + 1] + cnthi
         binct[i] <- binct[i] - cnthi
         binhi[i] <- binhi[i] - 1
         binlo[i + 1] <- binlo[i + 1] - 1
         changed = TRUE
      }
      else if (bestside == 0)
      {
         # move left
         binct[i - 1] <- binct[i - 1] + cntlo
         binct[i] <- binct[i] - cntlo
         binlo[i] <- binlo[i] + 1
         binhi[i - 1] <- binhi[i - 1] + 1
         changed = TRUE
      }
      
   }
   else
   {
      if (verbose) print("bins.move: no change")
   }
   
   names(binct) <- paste("[", xval[binlo], ", ", xval[binhi], "]", sep="")
   if (verbose) print(binct)
   return(list(binlo=binlo, binhi=binhi, binct=binct, xtbl=xtbl, xval=xval, changed=changed, 
               err=sqrt(bins.merr(binct, target.bins)), imax=bestindx, iside=bestside, gain=max(0, bestgain)))
}

#-------------------------------------------------------------------------------

#' \code{bins.split} - Split a bin into two bins optimally. 
#' 
#' @param force When \code{TRUE}, splits or merges bins regardless of whether the best gain is positive.
#' @param verbose When \code{TRUE}, prints resulting \code{binct}.
#' @export
#' @rdname bins.optimize
bins.split <- function(xval, xtbl, binlo, binhi, binct, target.bins, force=F, verbose=F)
{
   # try every break on every bin, pick the best
   nbins <- length(binlo)
   if (nbins != length(binhi) || nbins != length(binct)) stop("bins.split: lengths of binlo, binhi, binct do not match")
   err <- bins.merr(binct, target.bins)
   n2 <- (sum(binct) / target.bins)^2
   bestgain <- -Inf
   bestitem <- NA
   bestindx <- -1
   for (i in 1:nbins)
   {
      # binlo[i] = index of the smallest item in the bin
      # binhi[i] = index of the largest  item in the bin
      nitems <- binhi[i] - binlo[i] + 1
      if (nitems > 1)
      {  
         a <- 0
         b <- binct[i]
         for (k in binlo[i]:(binhi[i]-1))
         {
            a <- a + xtbl[k]
            b <- binct[i] - a
            gain   <- (2 * a * b + err - n2) / (nbins + 1)
            #print(paste(i, k, nitems, a, b, err, n2, gain))
            if (gain > bestgain)
            {
               bestgain <- gain
               bestitem <- k
               bestindx <- i
            }
         }
      }
   }
   
   # This works with a single bin too.
   changed = FALSE
   if (!is.na(bestitem) && (bestgain > 0 || force))
   {
      k <- bestitem
      a <- sum(xtbl[binlo[bestindx]:k])
      b <- binct[bestindx] - a
      binlo.new <- c(binlo[1:bestindx], k + 1)
      binhi.new <- c(k, binhi[bestindx:nbins])
      binct.new <- c(a, b)
      if (bestindx < nbins)
      {
         binlo.new <- c(binlo.new, binlo[(bestindx+1):nbins])
         binct.new <- c(binct.new, binct[(bestindx+1):nbins])
      }
      if (bestindx > 1)
      {
         binhi.new <- c(binhi[1:(bestindx-1)], binhi.new)
         binct.new <- c(binct[1:(bestindx-1)], binct.new)
      }
      binlo <- binlo.new
      binhi <- binhi.new
      binct <- binct.new      
      changed = TRUE
   }
   
   names(binct) <- paste("[", xval[binlo], ", ", xval[binhi], "]", sep="")
   if (verbose) print(binct)
   return(list(binlo=binlo, binhi=binhi, binct=binct, xtbl=xtbl, xval=xval, changed=changed,
               err=sqrt(bins.merr(binct, target.bins)), imax=bestindx, item=bestitem, gain=max(0, bestgain)))
}

#-------------------------------------------------------------------------------

#' \code{bins.merge} - Merges the two bins yielding the largest gain in error reduction. 
#' @export
#' @rdname bins.optimize
bins.merge <- function(xval, xtbl, binlo, binhi, binct, target.bins, force=F, verbose=F)
{
   # try every boundary, pick the best
   nbins <- length(binlo)
   if (nbins != length(binhi) || nbins != length(binct)) stop("bins.split: lengths of binlo, binhi, binct do not match")
   if (nbins < 2)
   {
      names(binct) <- paste("[", xval[binlo], ", ", xval[binhi], "]", sep="")
      if (verbose) print(binct)
      return(list(binlo=binlo, binhi=binhi, binct=binct, xtbl=xtbl, xval=xval, changed=FALSE,
                  err=sqrt(bins.merr(binct, target.bins)), imax=-1, gain=0))
   }
   err <- bins.merr(binct, target.bins)
   n2 <- (sum(binct) / target.bins)^2
   bestgain <- -Inf
   bestindx <- NA
   for (i in 1:(nbins - 1))
   {
      # take points from this bin if possible and place in the neighboring bin
      # binlo[i] = index of the smallest item in the bin
      # binhi[i] = index of the largest  item in the bin
      a <- binct[i]
      b <- binct[i + 1]
      gain <- -(2 * a * b + err - n2) / (nbins - 1)
      #print(paste(i, a, b, err, n2, gain))
      if (gain > bestgain)
      {
         bestgain <- gain
         bestindx <- i
      }
   }
   
   changed = FALSE
   if (!is.na(bestindx) && (bestgain > 0 || force))
   {
      aplusb <- binct[bestindx] + binct[bestindx + 1]
      binlo.new <- binlo[1:bestindx]
      binhi.new <- binhi[(bestindx + 1):nbins]
      binct.new <- aplusb
      if (bestindx < nbins - 1)
      {
         binlo.new <- c(binlo.new, binlo[(bestindx+2):nbins])
         binct.new <- c(binct.new, binct[(bestindx+2):nbins])
      }
      if (bestindx > 1)
      {
         binhi.new <- c(binhi[1:(bestindx-1)], binhi.new)
         binct.new <- c(binct[1:(bestindx-1)], binct.new)
      }
      binlo <- binlo.new
      binhi <- binhi.new
      binct <- binct.new      
      changed = TRUE
   }
   
   names(binct) <- paste("[", xval[binlo], ", ", xval[binhi], "]", sep="")
   if (verbose) print(binct)
   return(list(binlo=binlo, binhi=binhi, binct=binct, xtbl=xtbl, xval=xval, changed=changed,
               err=sqrt(bins.merr(binct, target.bins)), imax=bestindx, gain=max(0, bestgain)))
}

#-------------------------------------------------------------------------------

#' \code{bins.move.iter} - Apply \code{bins.move} until there's no change. Can only reduce the error.
#' @param lst List containing \code{xval, xtbl, binlo, binhi, binct}.
#' @rdname bins.optimize
bins.move.iter <- function(lst, target.bins, verbose=F)
{
   changed = TRUE
   while (changed) 
   { 
      lst <- bins.move(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, target.bins=target.bins, verbose=F)
      if (verbose) print(paste(lst$gain, lst$err))
      changed <- lst$changed
   } 
   return(lst)
}

#-------------------------------------------------------------------------------

#' \code{bins.split.iter}  Iterate to repeatedly apply \code{bins.split}. 
#' @param exact.groups If \code{FALSE}, run until either the target.bins is 
#' reached or there's no more splits or merges that reduce the error.
#' Otherwise (\code{TRUE}), run until the target.bins is reached, even if that 
#' increases the error.
#' @export
#' @rdname bins.optimize
bins.split.iter <- function(lst, target.bins, exact.groups=F, verbose=F)
{
   changed = TRUE
   # only split if the number of bins is less than desired
   while ( changed && length(lst$binct) < target.bins )  
   { 
      lst <- bins.split(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, target.bins=target.bins, force=exact.groups, verbose=verbose)
      if (verbose) print(paste(lst$gain, lst$err)) 
      changed <- lst$changed
   } 
   return(lst)
}

#-------------------------------------------------------------------------------

#' \code{bins.merge.iter} Iterate to repeatedly apply \code{bins.merge}. 
#' @export
#' @rdname bins.optimize
bins.merge.iter <- function(lst, target.bins, exact.groups=F, verbose=F)
{
   changed = TRUE
   # always merge if the number of bins exceeds the desired number
   while ( (!exact.groups && changed) || length(lst$binct) > target.bins ) 
   { 
      force = length(lst$binct) > target.bins
      lst <- bins.merge(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, target.bins=target.bins, force=force, verbose=verbose)
      if (verbose) print(paste(lst$gain, lst$err)) 
      changed <- lst$changed
   } 
   return(lst)
}

#-------------------------------------------------------------------------------
