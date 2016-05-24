# Rattle: A GUI for Data Mining in R
#
# BIN DATA
#
# Gnome R Data Miner: GNOME interface to R for Data Mining
#
# Time-stamp: <2014-09-05 21:27:32 gjw>
#
# Copyright (c) 2009-2014 Togaware Pty Ltd
#
# This files is part of Rattle.
#
# Rattle is free software: you can redistribute it and/or modify it
# under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 2 of the License, or
# (at your option) any later version.
#
# Rattle is distributed in the hope that it will be useful, but
# WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
# General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with Rattle. If not, see <http://www.gnu.org/licenses/>.
#
#-----------------------------------------------------------------------
#
# 070131 From Daniele Medri.
# 111025 Support wtd.quantiles suggested by Brenton R. Stone.

binning <- function (x, bins=4,
                     method=c("quantile", "wtd.quantile", "kmeans"),
                     labels=NULL, ordered=TRUE,
                     weights=NULL)
{
  # Set ordered to FALSE in Rattle since randomForests don't work when
  # the factor is ordered, for some reason (080406).
  
  # Best k for natural breaks

  varkmeans <- function (x, centers, iter.max=10, num.seeds=bins)
  {
    if (mode(x) == "numeric")
    {
      x <- data.frame(new.x=x)
    }
    KM <- kmeans(x=x, centers=centers, iter.max=iter.max)
    for (i in seq_len(num.seeds))
    {
      newKM <- kmeans(x=x, centers=centers, iter.max=iter.max)
      if (sum(newKM$withinss) < sum(KM$withinss))
      {
        KM <- newKM
      }
    }
    KM$tot.withinss <- sum(KM$withinss)
    xmean <- apply(x, 2, mean)
    centers <- rbind(KM$centers, xmean)
    bss1 <- as.matrix(dist(centers)^2)
    KM$betweenss <- sum(as.vector(bss1[nrow(bss1), ]) * c(KM$size, 0))
    return(KM)
  }

  method <- match.arg(method)
  if(is.factor(x)) stop(Rtxt("This variable is already a factor."))
  if (is.data.frame(x)) stop(Rtxt("An object of class data.frame is required."))
  if (length(x) < bins) stop(Rtxt("There are more bins than observations."))
  if (method == "wtd.quantile" &&
      ! packageIsAvailable("Hmisc", Rtxt("weighted quantile binning")))
    stop(Rtxt("wtd.quantile requires the Hmisc package."))
  
  # Binning

  x <- if (method == "quantile")
  {
    breaks <- c(quantile(x, probs = seq(0, 1, 1/bins), na.rm = TRUE, type=8))
    breaks <- unique(breaks)
    breaks[1] <- min(x, na.rm=TRUE)
    breaks[length(breaks)] <- max(x, na.rm=TRUE)
    # quantiles from quantile() can be non-unique, which cut() doesn't
    # like. This is handled above through unique(). The function
    # cut2() in Hmisc handles this situation gracefully and it could
    # be used, but it is not necessary.
    if(length(breaks) >= 2)
    {
      cut(x, breaks, include.lowest = TRUE, labels = labels)
    }
    else
    {
      cat(Rtxt("Warning: the variable is not considered.\n"))
      return(NULL)
    }
  }
  else if (method == "wtd.quantile")
  {
    breaks <- c(Hmisc::wtd.quantile(x, weights=weights, probs=seq(0, 1, 1/bins),
                                    na.rm=TRUE, type="quantile"))
    breaks <- unique(breaks)
    breaks[1] <- min(x, na.rm=TRUE)
    breaks[length(breaks)] <- max(x, na.rm=TRUE)
    # quantiles from quantile() can be non-unique, which cut() doesn't
    # like. This is handled above through unique(). The function
    # cut2() in Hmisc handles this situation gracefully and it could
    # be used, but it is not necessary.
    if(length(breaks) >= 2)
    {
      cut(x, breaks, include.lowest = TRUE, labels = labels)
    }
    else
    {
      cat(Rtxt("Warning: the variable is not considered.\n"))
      return(NULL)
    }
  }
  else if (method == "kmeans")
  {
    xx <- na.omit(x)
    maxbins <-nlevels(as.factor(xx))
    if(maxbins < bins)
    { 
      bins <-maxbins
    }
    breaks <- c(min(xx), tapply(xx, varkmeans(xx, bins)$cluster, max))
    if (length(unique(breaks)) >= 2)
    {
      cut(x, unique(breaks), include.lowest = TRUE, labels = labels)	
    }
    else
    {
      cat(Rtxt("Warning: the variable is not considered.\n"))
      return(NULL)	
    }
  }

  if(ordered == TRUE)
    result <- ordered(factor(x))
  else
    result <- factor(x)

  attr(result, "breaks") <- breaks
  return(result)
}
