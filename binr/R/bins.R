#-------------------------------------------------------------------------------
#
# Package binr 
#
# Implementation. 
# 
# Sergei Izrailev, 2011-2015
#-------------------------------------------------------------------------------
# Copyright 2011-2014 Collective, Inc.
# Copyright 2015 Jabiru Ventures LLC
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


#' Cut Numeric Values Into Evenly Distributed Groups (Bins)
#' 
#' \code{bins} - Cuts points in vector x into evenly distributed groups (bins). 
#' \code{bins} takes 3 separate approaches to generating the cuts, picks the one 
#' resulting in the least mean square deviation from the ideal cut -
#' \code{length(x) / target.bins} points in each bin - and then  merges small bins 
#' unless excat.groups is \code{TRUE}
#' The 3 approaches are: 
#' \enumerate{
#' \item{Use quantiles, and increase the number of even cuts up to max.breaks until the 
#'    number of groups reaches the desired number. See \code{\link{bins.quantiles}}. } 
#' \item{Start with a single bin with all the data in it and perform bin splits until 
#'    either the desired number of bins is reached or there's no reduction in error
#'    (the latter is ignored if \code{exact.groups} is \code{TRUE}). See \code{\link{bins.split}}. }
#' \item{Start with \code{length(table(x))} bins, each containing exactly one distinct value 
#'    and merge bins until the desired number of bins is reached. If \code{exact.groups} is
#'    \code{FALSE}, continue merging until there's no further reduction in error. 
#'    See \code{\link{bins.merge}}. }
#' }
#' For each of these approaches, apply redistribution of points among existing bins
#' until there's no further decrease in error. See \code{bins.move}.
#'
#' The gains are computed using incremental analytical expresions derived for moving 
#' a value from one bin to the next, splitting a bin into two or merging two bins.
#' @param x Vector of numbers
#' @param target.bins Number of groups desired; this is also the max number of groups.
#' @param max.breaks   Used for initial cut. If \code{exact.groups} is \code{FALSE}, bins are merged 
#'                     until there's no bins with fewer than \code{length(x) / max.breaks} points.
#'                     In \code{bins}, one of \code{max.breaks} and \code{minpts} must be supplied.
#' @param exact.groups if TRUE, the result will have exactly the number of target.bins;
#'                     if FALSE, the result may contain fewer than target.bins bins
#' @param verbose      Indicates verbose output.
#' @param errthresh    If the error is below the provided value, stops after the first rough estimate of the bins.
#' @param minpts       Minimum number of points in a bin. 
#'                     In \code{bins}, one of \code{max.breaks} and \code{minpts} must be supplied.
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
#' 
#' @examples
#' \dontrun{
#'    # Seriously skewed x:
#'    x <- floor(exp(rnorm(200000 * 1.3)))
#'    cuts <- bins(x, target.bins = 10, minpts = 2000)
#'    cuts$breaks <- bins.getvals(cuts)
#'    cuts$binct
#'    #   [0, 0]    [1, 1]    [2, 2]    [3, 3]    [4, 4]    [5, 5]    [6, 7]   [8, 10] 
#'    # 129868     66611     28039     13757      7595      4550      4623      2791 
#'    #   [11, 199] 
#'    # 2166 
#' 
#'    # Centered x:
#'    x <- rep(c(1:10,20,31:40), c(rep(1, 10), 100, rep(1,10)))
#'    cuts <- bins(x, target.bins = 3, minpts = 10)
#'    cuts$binct
#'    # [1, 10] [20, 20] [31, 40] 
#'    #      10      100       10 
#' }
#' @seealso \code{\link{binr}}, \code{\link{bins.greedy}}, \code{\link{bins.quantiles}} \code{\link{bins.optimize}}
#' @export
#' @rdname bins
bins <- function(x, target.bins, max.breaks = NA, exact.groups=F, verbose=F, errthresh = 0.1, minpts = NA)
{
   if (length(x) < target.bins) stop(paste("bins: number of desired groups (", target.bins, ") is greater than the number of points (", length(x), ")"))
   all <- vector(3, mode="list")
   names(all) <- c("quantile", "split", "merge")
   if (is.na(minpts)) 
   { 
      if (is.na(max.breaks) || max.breaks < 1)
      {
         stop(paste("bins: missing or invalid max.breaks = ", max.breaks, " and minpts = ", minpts))
      }
      minpts <- floor(length(x) / max.breaks)
   } else if (is.na(max.breaks)) 
   {
      if (minpts < 1) stop(paste("bins: invalid minpts = ", minpts))
      max.breaks <- ceiling(length(x) / minpts)
   } else
   {
      stop("bins: only one of max.breaks and minpts can be specified")
   }
   
   all$quantile <- bins.quantiles(x, target.bins, max.breaks)
   if (all$quantile$err / (length(x) / target.bins) < errthresh) return(all$quantile) # good enough
   
   xval <- all$quantile$xval
   xtbl <- all$quantile$xtbl
   
   all$split <- list(binlo = 1, binhi = length(xval), binct = sum(xtbl), xval = xval, xtbl = xtbl)
   #all$merge <- list(binlo = 1:length(xtbl), binhi = 1:length(xtbl), binct = xtbl, xval = xval, xtbl = xtbl)   
   all$merge <- bins.quantiles(x, max.breaks * 2, max.breaks * 2) # start with a fairly fine breakup 
   
   all$quantile <- bins.move.iter(all$quantile, target.bins, verbose=verbose)
   if (exact.groups)
   {
      if (length(all$quantile$binct) < exact.groups)
      {
         all$quantile <- bins.split.iter(all$quantile, target.bins, verbose=verbose, exact.groups=T)         
         all$quantile <- bins.move.iter(all$quantile, target.bins, verbose=verbose)
      }
      else if (length(all$quantile$binct) > exact.groups)
      {
         all$quantile <- bins.merge.iter(all$quantile, target.bins, verbose=verbose, exact.groups=T)         
         all$quantile <- bins.move.iter(all$quantile, target.bins, verbose=verbose)
      }
   }
   
   # don't attempt to split or merge if the number of values is too big
   if (length(xval) < 10000)
   {
      # split
      all$split <- bins.split.iter(all$split, target.bins, verbose=verbose, exact.groups=exact.groups)
      all$merge <- bins.merge.iter(all$merge, target.bins, verbose=verbose, exact.groups=exact.groups)      
   
      all <- lapply(all, function(lst) bins.move.iter(lst, target.bins, verbose=verbose))
      errs <- unlist(lapply(all, function(lst) sqrt(bins.merr(lst$binct, target.bins))))
      names(errs) <- names(all)
      if (verbose) print(errs)
      imin <- which(errs == min(errs))[1]
      
      lst <- all[[imin]]
   }
   else
   {
      lst <- bins.move.iter(all$quantile, target.bins, verbose=verbose)
   }

   if (!exact.groups)
   {
      # merge small bins
      while (any(lst$binct < minpts))
      {
         lst <- bins.merge(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, target.bins=target.bins, force=T, verbose=verbose)
         lst <- bins.move.iter(lst, target.bins, verbose=verbose)
         if (verbose) print(paste(lst$gain, lst$err))          
      }
   }
   
   return(lst)
}


#-------------------------------------------------------------------------------

#' \code{bins.getvals} - Extracts cut points from the object retured by \code{bins}. 
#' The cut points are always between the values in \code{x} and weighed such that 
#' the cut point splits the area under the line from (lo, n1) to (hi, n2) in half.
#' @param lst The list returned by the \code{bins} function.
#' @param minpt The value replacing the lower bound of the cut points.
#' @param maxpt The value replacing the upper bound of the cut points.
#' @return \code{bins.getvals} returns a vector of cut points extracted from the 
#'         \code{lst} object. 
#' @export
#' @rdname bins
bins.getvals <- function(lst, minpt = -Inf, maxpt = Inf)
{
   # finds the point that splits the area under the line from (lo, n1) to (hi, n2) in half.
   cutpt <- function(lo, hi, n1, n2)
   {
      if (n1 == n2) return((hi + lo) / 2)
      a <- hi - lo
      b <- 2 * a * n1 / (n2 - n1)
      c = a * a * (n1 + n2)  / (n2 - n1) / 2
      d <- (-b + sign(n2 - n1) * sqrt(b * b + 4 * c)) / 2
      return(lo + d)
   }
   
   nbins <- length(lst$binct)
   if (nbins == 0) stop("bins.getvals: zero bins")
   res <- vector(nbins + 1, mode="double")
   res[1] <- minpt
   res[nbins + 1] <- maxpt
   
   if (nbins > 1)
   {
      for (i in 2:nbins)
      {
         res[i] <- cutpt(lst$xval[lst$binhi[i - 1]], lst$xval[lst$binlo[i]], lst$binct[i - 1], lst$binct[i])
      }
      names(res) <- c(paste("[", res[1:(nbins-1)], ", ", res[2:nbins], ")", sep=""), paste("[", res[nbins], ", ", res[nbins + 1], "]", sep=""), "LAST")     
   } else {
      # nbins = 1
      names(res) <- c(paste("[", res[1], ", ", res[2], "]", sep=""), "LAST")
   }
   
   attr(res, "binlo") <- lst$xval[lst$binlo]
   attr(res, "binhi") <- lst$xval[lst$binhi]
   attr(res, "binct") <- lst$binct
   return(res)
}


#-------------------------------------------------------------------------------
 
#' \code{bins.merr} - Partitioning the data into bins using splitting, merging 
#' and moving optimizes this error function, which is the mean squared error 
#' of point counts in the bins relative to the optimal number of points per bin.
#' @param binct The number of points falling into the bins.
#' @export
#' @rdname bins
bins.merr <- function(binct, target.bins) { mean((binct - sum(binct) / target.bins)^2) }
   
#-------------------------------------------------------------------------------

if (0)
{
   xnorm <- rnorm(200000)
   xlog <- floor(exp(xnorm * 1.3))
   xlog2 <- c(xlog, -xlog - 1)
   x <- xlog2
   
   xtbl <- table(x)
   xval <- sort(unique(x))
   nvals <- length(xval)
   nbins <- 10
   
   binsz <- floor(length(x) / nbins)
   minpts <- 1000
   thresh <- 0.8
   target.bins = 10
   max.breaks = 20
   
   binlo <- 1:10
   binhi <- c(1:9,length(xval))
   binct <- c(xtbl[1:9],sum(xtbl[10:length(xtbl)]))
   source("~/cpp/cmrutils/R/cmcut.R")
   lst <- adjust.cells(xval, xtbl, binlo, binhi, binct)
   bc <- rep(0,length(lst$binct)); while(!all(lst$binct == bc)) { bc <- lst$binct; lst <- adjust.cells(xval, xtbl, lst$binlo, lst$binhi, lst$binct) }
   lst <- split.cells(xval, xtbl, lst$binlo, lst$binhi, lst$binct); lst
   prevlen <- 0; while(length(lst$binct) != prevlen) { prevlen <- length(lst$binct); lst <- split.cells(xval, xtbl, lst$binlo, lst$binhi, lst$binct) }
 
   lst <- bins.quantiles(x, 10, 20)
   preverr <- Inf; err <- sqrt(bins.merr(lst$binct, target.bins)); while (preverr > err) { preverr <- err; lst <- bins.move(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct); err <- sqrt(bins.merr(lst$binct, target.bins)); print(paste(lst$gain, err)); }
   
   lst$binlo <- 1:10
   lst$binhi <- c(1:9,length(xval))
   lst$binct <- c(xtbl[1:9],sum(xtbl[10:length(xtbl)]))
   
   lst$binlo <- 1; lst$binhi <- length(xval); lst$binct <- sum(xtbl)
   preverr <- Inf; err <- sqrt(bins.merr(lst$binct, target.bins)); while (preverr > err) { preverr <- err; lst <- bins.split(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, target.bins); err <- sqrt(bins.merr(lst$binct, target.bins)); print(paste(lst$gain, err)); }
   preverr <- Inf; err <- sqrt(bins.merr(lst$binct, target.bins)); while (preverr > err) { preverr <- err; lst <- bins.move(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct); err <- sqrt(bins.merr(lst$binct, target.bins)); print(paste(lst$gain, err)); }  
   preverr <- Inf; err <- sqrt(bins.merr(lst$binct, target.bins)); while (preverr > err) { preverr <- err; lst <- bins.merge(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, target.bins); err <- sqrt(bins.merr(lst$binct, target.bins)); print(paste(lst$gain, err)); }
   preverr <- Inf; err <- sqrt(bins.merr(lst$binct, target.bins)); while (preverr > err) { preverr <- err; lst <- bins.move(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct); err <- sqrt(bins.merr(lst$binct, target.bins)); print(paste(lst$gain, err)); }

   
   lst$binlo <- 1; lst$binhi <- length(xval); lst$binct <- sum(xtbl)
   preverr <- Inf; err <- sqrt(bins.merr(lst$binct, target.bins)); while (preverr > err) { preverr <- err; lst <- bins.split(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, 10); err <- sqrt(bins.merr(lst$binct, target.bins)); print(paste(lst$gain, err)); }
   while(length(lst$binct) < 10) { lst <- bins.split(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, target.bins, force=T) }
   preverr <- Inf; err <- sqrt(bins.merr(lst$binct, target.bins)); while (preverr > err) { preverr <- err; lst <- bins.split(lst$xval, lst$xtbl, lst$binlo, lst$binhi, lst$binct, 10); err <- sqrt(bins.merr(lst$binct, target.bins)); print(paste(lst$gain, err)); }
 }
 

#-------------------------------------------------------------------------------
