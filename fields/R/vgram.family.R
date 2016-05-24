# fields  is a package for analysis of spatial data written for
# the R software environment .
# Copyright (C) 2016
# University Corporation for Atmospheric Research (UCAR)
# Contact: Douglas Nychka, nychka@ucar.edu,
# National Center for Atmospheric Research, PO Box 3000, Boulder, CO 80307-3000
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with the R software environment if not, write to the Free Software
# Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301  USA
# or see http://www.r-project.org/Licenses/GPL-2    
"vgram" <- function(loc, y, id = NULL, d = NULL, lon.lat = FALSE, 
                    dmax = NULL, N = NULL, breaks = NULL, 
                    type=c("variogram", "covariogram", "correlogram")) {
  
  type=match.arg(type)
  
  # coerce to matrix
  y <- cbind(y)
  # if nearest neighbor indices are missing create all possible pairs.
  if (is.null(id)) {
    n <- nrow(loc)
    is = rep(1:n, n)
    js = rep(1:n, rep(n, n))
    ind <- is > js
    id <- cbind(is, js)[ind, ]
  }
  
  # if distances are missing calculate these
  if (is.null(d)) {
    loc <- as.matrix(loc)
    if (lon.lat) {
      d <- rdist.earth.vec(loc[id[,1],], loc[id[,2],]) #we want result in miles, not meters
    }
    else {
      d <- rdist.vec(loc[id[,1],], loc[id[,2],])
    }
  }
  
  # normalize columns to create correlogram, if necessary
  #
  if(type == "correlogram") {
    sigma = apply(y, 2, sd, na.rm=TRUE)
    y = sweep(y, 2, (1/sigma), FUN="*")
  }
  
  # center the columns by their mean and get row means if y is a matrix
  #
  colMeans <- apply(y, 2, mean, na.rm=TRUE)
  yCntr = sweep(y, 2, colMeans) 
  y1Cntr = yCntr[id[,1],]
  y2Cntr = yCntr[id[,2],]
  if(type == "variogram") {
    vg <- 0.5 * rowMeans(cbind((y1Cntr - y2Cntr)^2), 
                         na.rm = TRUE)
  }
  else {
    vg <- rowMeans(cbind(y1Cntr * y2Cntr), 
                   na.rm = TRUE)
  }
  #
  #information for returned object
  #
  call <- match.call()
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  ind <- d <= dmax & !is.na(vg)
  
  ## add a binned  variogram if breaks are supplied
  out <- list(d = d[ind], vgram = vg[ind], call = call, type=type)
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("vgram", class(out))
  out
}

#calculating cross-covariogram and cross-correlogram (cross-covariance and 
#cross-correlation)
crossCoVGram = function(loc1, loc2, y1, y2, id = NULL, d = NULL, lon.lat = FALSE, 
                        dmax = NULL, N = NULL, breaks = NULL, 
                        type=c("cross-covariogram", "cross-correlogram")) {
  
  type=match.arg(type)
  
  # coerce to matrix
  y1 <- cbind(y1)
  y2 <- cbind(y2)
  
  # if nearest neighbor indices are missing create all possible pairs.
  if (is.null(id)) {
    n1 <- nrow(loc1)
    n2 <- nrow(loc2)
    id <- cbind(rep(1:n1, n2), rep(1:n2, rep(n1, n2)))
  }
  
  # if distances are missing calculate these
  if (is.null(d)) {
    loc1 <- as.matrix(loc1)
    loc2 <- as.matrix(loc2)
    if (lon.lat) {
      d <- rdist.earth.vec(loc1[id[,1],], loc2[id[,2],]) #we want result in miles, not meters
    }
    else {
      d <- rdist.vec(loc1[id[,1],], loc2[id[,2],])
    }
  }
  #
  # calculating covariogram will center the columns by their mean and get row means if y is a matrix
  #
  colMeans1 <- apply(y1, 2, mean, na.rm=TRUE)
  colMeans2 <- apply(y2, 2, mean, na.rm=TRUE)
  y1Cntr = sweep(data.matrix(y1), 2, colMeans1)  # subtract the column means
  y2Cntr = sweep(data.matrix(y2), 2, colMeans2)  # subtract the column means
  #
  # normalize to create cross-correlogram, if necessary
  #
  if(type == "cross-correlogram") {
    sigma1 = apply(y1Cntr, 2, sd, na.rm=TRUE)
    sigma2 = apply(y2Cntr, 2, sd, na.rm=TRUE)
    y1Cntr = sweep(y1Cntr, 2, 1/sigma1, FUN="*")
    y2Cntr = sweep(y2Cntr, 2, 1/sigma2, FUN="*")
  }
  #
  # calculate covariance for the given points
  #
  y1Cntr = y1Cntr[id[,1],]
  y2Cntr = y2Cntr[id[,2],]
  vg <- rowMeans(cbind(y1Cntr*y2Cntr), na.rm = TRUE)
  #
  #information for returned object
  #
  call <- match.call()
  if (is.null(dmax)) {
    dmax <- max(d)
  }
  od <- order(d)
  d <- d[od]
  vg <- vg[od]
  ind <- d <= dmax & !is.na(vg)
  ## add a binned  variogram if breaks are supplied
  out <- list(d = d[ind], vgram = vg[ind], call = call, type=type)
  if (!is.null(breaks) | !is.null(N)) {
    out <- c(out, stats.bin(d[ind], vg[ind], N = N, breaks = breaks))
  }
  class(out) = c("vgram", class(out))
  out
}

#plot only the line of the empirical variogram, where the y coordinates of the line are 
#at the means of the bins
plot.vgram = function(x, N=10, breaks = pretty(x$d, N, eps.correct = 1), add=FALSE, ...) {
  otherArgs = list(...)
  type=x$type
  
  #set y axis label if not set by user
  if(is.null(otherArgs$ylab)) {
    if(type=="variogram")
      otherArgs$ylab = "sqrt(Variance)"
    else if(type == "covariogram" || type=="cross-covariogram")
      otherArgs$ylab = "Covariance"
    else if(type == "correlogram" || type=="cross-correlogram")
      otherArgs$ylab = "Correlation"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  
  #set x axis label if not set by user
  if(is.null(otherArgs$xlab))
    otherArgs$xlab = "Distance"
  
  #set plot title if not set by user
  if(is.null(otherArgs$main)) {
    if(type=="variogram")
      otherArgs$main = "Empirical Variogram"
    else if(type=="covariogram")
      otherArgs$main = "Empirical Covariogram"
    else if(type=="correlogram")
      otherArgs$main = "Empirical Correlogram"
    else if(type=="cross-covariogram")
      otherArgs$main = "Empirical Cross-Covariogram"
    else if(type=="cross-correlogram")
      otherArgs$main = "Empirical Cross-Correlogram"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  
  #set ylim for correlogram if not set by user
  if(is.null(otherArgs$ylim)) {
    if(type == "correlogram" || type=="cross-correlogram")
      otherArgs$ylim = c(-1, 1)
  }
  
  #set line type if not set by user
  if(is.null(otherArgs$type))
    otherArgs$type = "o"
  
  #get bin data
  dat = getVGMean(x, breaks=breaks)
  
  #get bin centers versus bin means
  centers = dat$centers
  ys = dat$ys
  
  #remove NAs
  notNas = !is.na(ys)
  centers = centers[notNas]
  ys = ys[notNas]
  
  #plot
  if(!add)
    do.call("plot", c(list(centers, ys), otherArgs))
  else
    do.call("lines", c(list(centers, ys), otherArgs))
}

"boxplotVGram" = function(x, N=10, breaks = pretty(x$d, N, eps.correct = 1), plot=TRUE, 
                          plot.args=NULL, ...) {
  dists = x$d
  type=x$type
  if(type == "variogram")
    y = sqrt(x$vgram)
  else
    y = x$vgram
  otherArgs = list(...)
  
  #set y axis label if not set by user
  if(is.null(otherArgs$ylab)) {
    if(type=="variogram")
      otherArgs$ylab = "sqrt(Variance)"
    else if(type == "covariogram" || type=="cross-covariogram")
      otherArgs$ylab = "Covariance"
    else if(type == "correlogram" || type=="cross-correlogram")
      otherArgs$ylab = "Correlation"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  
  #set x axis label if not set by user
  if(is.null(otherArgs$xlab))
    otherArgs$xlab = "Distance"
  
  #set plot title if not set by user
  if(is.null(otherArgs$main)) {
    if(type=="variogram")
      otherArgs$main = "Empirical Variogram"
    else if(type=="covariogram")
      otherArgs$main = "Empirical Covariogram"
    else if(type=="correlogram")
      otherArgs$main = "Empirical Correlogram"
    else if(type=="cross-covariogram")
      otherArgs$main = "Empirical Cross-Covariogram"
    else if(type=="cross-correlogram")
      otherArgs$main = "Empirical Cross-Correlogram"
    else
      stop("vgram 'type' argument must be either 'variogram', 'covariogram', 'correlogram', 'cross-covariogram', or 'cross-correlogram'")
  }
  
  #set ylim for correlogram if not set by user
  if(is.null(otherArgs$ylim)) {
    if(type == "correlogram" || type=="cross-correlogram")
      otherArgs$ylim = c(-1, 1)
  }
  
  #make boxplot
  bplot = do.call("bplot.xy", c(list(x=dists, y=y, N=N, breaks=breaks, plot=plot), otherArgs))
  
  #return bplot.xy statistics if plot==FALSE
  if(!plot)
    return(bplot)
  
  #plot bin means with plot parameters given in plot.args (with defaults to look pretty)
  plot.args$x=x
  plot.args$add=TRUE
  plot.args$breaks=breaks
  if(is.null(plot.args$col))
    plot.args$col = "red"
  if(is.null(plot.args$type))
    plot.args$type = "p"
  do.call("plot.vgram", plot.args)
}

# Returns the variogram bin centers and means
getVGMean = function(x, N = 10, breaks = pretty(x$d, N, eps.correct = 1)) 
{
  # Can calculate mean or other statistical functions of the values in the bins
  VGstat = function(VG, minD=-Inf, maxD=Inf, statFun="mean", ...) {
    ind = (VG$d > minD) & (VG$d < maxD)
    do.call(statFun, c(list(VG$vgram[ind]), list(...)))
  }
  
  #helper function to get mean from any single bin
  meansFromBreak = function(breakBounds = c(-Inf, Inf)) {
    VGstat(x, minD=breakBounds[1], maxD=breakBounds[2], na.rm=TRUE)
  }
  
  #apply helper function to all bins
  lowBreaks = breaks
  highBreaks = c(breaks[2:length(breaks)], Inf)
  breakBounds = cbind(lowBreaks, highBreaks)
  centers = apply(breakBounds, 1, mean, na.rm=TRUE)
  ys = apply(breakBounds, 1, meansFromBreak)
  
  #take square root if variogram
  if(x$type == "variogram")
    ys=sqrt(ys)
  
  return(list(centers=centers, ys=ys, type=x$type))
}