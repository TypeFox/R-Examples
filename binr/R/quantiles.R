#-------------------------------------------------------------------------------
#
# Package binr 
#
# Implementation of an algorithm based on quantiles. 
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


#' Cuts the data set x into roughly equal groups using quantiles. 
#' 
#' Because the number of unique values may be smaller than target.bins,
#' the function gradually increases the number of quantiles up to max.breaks 
#' or until the target.bins number of bins is reached.
#' @param x A numeric vector to be cut in bins.
#' @param target.bins Target number of bins, which may not be reached if the 
#'        number of unique values is smaller than the specified value.
#' @param max.breaks Maximum number of quantiles; must be at least as 
#'        large as \code{target.bins}. 
#' @param verbose Indicates verbose output.
#' @name bins.quantiles
#' @title Quantile-based binning
#' @seealso \code{\link{binr}}, \code{\link{bins}}, \code{\link{bins.greedy}}, \code{\link{bins.optimize}}
#' @export
bins.quantiles <- function(x, target.bins, max.breaks, verbose = FALSE) 
{
   n.groups <- 1
   nbreak <- target.bins - 1 
   if (max.breaks < target.bins) stop("bins.quantiles: max.breaks must be at least as large as target.bins")
   xtbl <- table(x)
   xval <- sort(unique(x))
   while (n.groups <= target.bins && nbreak <= max.breaks) 
   {
      nbreak <- nbreak + 1
      qq <- unique(quantile(x, probs = 0:nbreak / nbreak))
      idx <- unique(findInterval(qq, xval))
      n.groups <- length(idx)
   }
   idx <- unique(findInterval(qq, xval))
   nbins <- length(idx) - 1
   binlo <- vector(nbins, mode="integer")
   binhi <- vector(nbins, mode="integer")
   binct <- vector(nbins, mode="integer")
   prev <- 0
   for (i in 1:nbins)
   {
      binlo[i] <- prev + 1
      binhi[i] <- idx[i + 1]
      binct[i] <- sum(xtbl[binlo[i]:binhi[i]])
      prev <- idx[i + 1]
   }
   names(binct) <- paste("[", xval[binlo], ", ", xval[binhi], "]", sep="")
   if (verbose) print(binct)
   return(list(binlo=binlo, binhi=binhi, binct=binct, xtbl=xtbl, xval=xval, err=sqrt(bins.merr(binct, target.bins))))
}

#-------------------------------------------------------------------------------
