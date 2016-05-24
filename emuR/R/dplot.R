##' A function to plot one or more columns of EMU-trackdata as a function of
##' time
##' 
##' A general purpose routine for plotting EMU-trackdata on a single plot.
##' Tracks can be aligned at an arbitrary position, length normalised or
##' averaged. The plots can be colour-coded for different category types.
##' 
##' 
##' @aliases dplot dplot.norm dplot.time
##' @param x An EMU-trackdata object
##' @param labs A label vector with one element for each row in 'dataset'
##' @param offset Either: A single numeric vector between 0 and 1. 0 and 1
##' denote synchronize the trackdata at their temporal onsets and offsets
##' respectively; 0.5 denotes synchronization at the temporal midpoint, etc. Or
##' a numeric vector of the same length as x specifying the synchronisation
##' point per segment
##' @param prop A single element character vector specifying whether the tracks
##' should be aligned proportionally or relative to millisecond times. Defaults
##' to proportional alignment
##' @param average If TRUE, the data for each unique label in 'labs' is
##' averaged
##' @param xlim A vector of two numeric values specifying the x-axis range
##' @param ylim A vector of two numeric values specifying the y-axis range
##' @param lty A single element logical vector. Defaults to F.  If TRUE, plot
##' each label type in a different linetype
##' @param normalise If TRUE, the data for each segment is linearly time
##' normalised so that all observations have the same length. The number of
##' points used in the linear time normalisation is control by the argument n.
##' @param colour A single element logical vector. Defaults to T to plot each
##' label type in a different colour
##' @param lwd A code passed to the lwd argument in plotting functions. 'lwd'
##' can be either a single element numeric vector, or its length must be equal
##' to the number of unique types in labs. For example, if lwd=3 and if labs =
##' c("a", "b", "a", "c"), then the output is c(3, 3, 3, 3). Alternatively, if
##' lwd = c(2,3,1), then the output is c(2, 3, 2, 1) for the same example. The
##' default is NULL in which case all lines are drawn with lwd=1
##' @param pch A code passed to the pch argument in plotting functions.
##' Functions in the same way as lwd above
##' @param legend Either a character vector to plot the legend. Possible values
##' are: "bottomright"', '"bottom"', '"bottomleft"', '"left"', '"topleft"',
##' '"top"', '"topright"', '"right"' and '"center"'. This places the legend on
##' the inside of the plot frame at the given location. Partial argument
##' matching is used. Or a logical vector: legend = FALSE suppresses legend
##' plotting. legend = TRUE plots it at the default, legend = "topright"
##' @param axes A single element logical vector. Defaults to T to plot the axes
##' @param type The default line type. Default to "l" for a line plot
##' @param n A single element numeric vector. Only used if normalise=T. The
##' number of data points used to linearly time normalise each track
##' @param ... graphical options \link{par}
##' @return NULL
##' @author Jonathan Harrington
##' @seealso \code{\link{dcut}} \code{\link{get_trackdata}}
##' @keywords dplot
##' @examples
##' 
##' 
##'    # Plot of column 1 (which happens to be the 1st formant) of an EMU-trackdata object
##'    dplot(dip.fdat[,1])
##' 	
##' 
##'    # As above but only observations 1 to 5
##'    dplot(dip.fdat[1:5,1])
##' 	
##' 
##'    #  column 2 (which happens to be of the second formant) and colour-coded
##'    # for each label-type
##'    dplot(dip.fdat[,2], dip.l)
##' 	
##' 
##'    # put the legend bottom left
##'    dplot(dip.fdat[,2], dip.l, legend="bottomleft")
##' 	
##' 
##'    # as above with no legend and averaged per category
##'    dplot(dip.fdat[,2], dip.l, legend=FALSE, average=TRUE)
##' 	
##' 
##'    # both formants averaged
##'    dplot(dip.fdat[,1:2], dip.l, average=TRUE)
##' 	
##' 
##'    # F2 only with linear-time normalisation
##'    dplot(dip.fdat[,2], dip.l, norm=TRUE)
##' 	
##' 
##'    # linear time-normalisation, both formants and averaged
##'    dplot(dip.fdat[,1:2], dip.l, norm=TRUE, average=TRUE)
##' 	
##' 
##'    # synchronise at the temporal midpoint before averaging, F2 only
##'    dplot(dip.fdat[,2], dip.l, offset=0.5, average=TRUE)
##' 	
##' 
##'    # synchronise 60 ms before the diphthong offset
##'    dplot(dip.fdat[,2], dip.l, offset=dip.fdat$ftime[,2]-60, prop=FALSE)
##' 	
##' 
##'    # as above averaged, no colour with linetype, 
##' # different plot symbols double line thickness in the range between +- 20 ms
##'    dplot(dip.fdat[,2], dip.l, offset=dip.fdat$ftime[,2]-60, prop=FALSE,
##'    average=TRUE, colour=FALSE, lty=TRUE, pch=1:3, lwd =2, type="b", xlim=c(-20, 20))
##' 
##' 
##' 
##' 
##' @export dplot
`dplot` <- function (x, labs = NULL, offset = 0, prop=TRUE, 
                     average = FALSE, xlim = NULL, ylim = NULL,  
                     lty = FALSE, normalise = FALSE, colour = TRUE, 
                     lwd = NULL, pch=NULL, legend = "topright", 
                     axes = TRUE, type="l", n = 20, ...) 
{
  
  
  if(prop)
  {
    if(length(offset) != 1)
      stop("Specify only one offset time when prop=T")
    else if ((offset < 0) | (offset > 1)) 
      stop("offset must be between 0 and 1 when prop=T")
    
    
  }
  
  else
    if(nrow(x) != length(offset))
      stop("nrow(x) and length(offset) must be the same when prop=F")
  
  pout <- NULL
  if (is.matrix(x$data)) {
    pout <- as.list(NULL)
    pout$data <- as.list(NULL)
    mat <- NULL
    if (is.null(ylim)) 
      ylim <- range(x$data)
    numcols <- ncol(x$data)
    
    for (j in 1:ncol(x$data)) {
      mat <- x
      mat$data <- mat$data[, j]
      if (!normalise) 
        vals <- dplot.time(mat, labs = labs, offset = offset, 
                           prop=prop, average = average,  xlim = xlim, 
                           ylim = ylim,  lty = lty, 
                           colour = colour, legend = legend, lwd = lwd, pch=pch, type=type)
      else vals <- dplot.norm(mat, labs = labs, average = average, 
                              xlim = xlim, ylim = ylim,  lty = lty, 
                              colour = colour, legend = legend, lwd = lwd, pch=pch, type=type,  n = n)
      graphics::par(new = TRUE)
      pout$data[[j]] <- vals$data
      if (j == ncol(x$data)) {
        pout$time <- vals$time
        pout$labs <- vals$labs
      }
    }
  }
  else {
    if (!normalise) 
      pout <- dplot.time(x, labs = labs, offset = offset, 
                         prop=prop, average = average,  xlim = xlim, ylim = ylim, 
                         lty = lty, colour = colour, 
                         lwd = lwd, pch=pch, type=type, legend = legend)
    else pout <- dplot.norm(x, labs = labs, average = average, xlim = xlim, 
                            ylim = ylim,  lty = lty, colour = colour, 
                            lwd = lwd, pch=pch, type=type, legend = legend,  n = n)
  }
  
  graphics::par(new = FALSE)
  
  invisible(pout)
  graphics::title(...)
  graphics::box(...)
  if (axes) {
    graphics::axis(side = 1)
    graphics::axis(side = 2)
  }
}



##' @export
`dplot.time` <- function (x, labs = NULL, offset = 0, prop=TRUE, 
                          average = FALSE,  xlim = NULL, ylim = NULL, 
                          lty = FALSE, colour = TRUE, lwd = NULL, 
                          pch=NULL, legend = "topright",  type="l", ...) 
{
  ovec <- as.list(NULL)
  samrate <- 1000/((x$ftime[1, 2] - x$ftime[1, 
                                            1])/(x$index[1, 2] - x$index[1, 1]))
  if (is.null(labs)) 
    labs <- rep(1, nrow(x$index))
  col.lty <- mu.colour(labs, colour, lty, lwd, pch)
  colour <- col.lty$colour
  lty <- col.lty$linetype
  lwd <- col.lty$lwd
  pch <- col.lty$pch
  if (prop) 
    ref.time <- x$ftime[, 1] + ((x$ftime[, 2] - 
                                   x$ftime[, 1]) * offset)
  else ref.time <- offset
  maxlen <- 2 * (max(x$index[, 2] - x$index[, 1] + 
                       1))
  pointval <- round(maxlen/2)
  mat.na <- matrix(NA, nrow(x$index), maxlen)
  for (j in 1:nrow(x$index)) {
    left <- x$index[j, 1]
    right <- x$index[j, 2]
    length.index <- right - left + 1
    times <- x$ftime[j, ]
    refn <- ref.time[j]
    inval <- closest(seq(times[1], times[2], length = length.index), 
                     refn)
    inval <- inval[1]
    left.na <- pointval - inval + 1
    right.na <- left.na + length.index - 1
    mat.na[j, left.na:right.na] <- x$data[left:right]
  }
  z <- apply(mat.na, 2, mean, na.rm = TRUE)
  natemp <- is.na(z)
  nums <- c(1:length(natemp))
  nonums <- nums[!natemp]
  interval <- 1000/samrate
  if (is.null(xlim)) 
    xlim <- c(nonums[1], nonums[length(nonums)])
  else xlim <- c(pointval + xlim[1]/interval, pointval + xlim[2]/interval)
  time1 <- (1 - pointval) * interval
  time2 <- (ncol(mat.na) - pointval) * interval
  xtime <- seq(time1, time2, length = ncol(mat.na))
  xtimelim <- (xlim - pointval) * interval
  if (is.null(ylim)) 
    ylim <- range(mat.na, na.rm = TRUE)
  if (!average) {
    for (j in 1:nrow(mat.na)) {
      graphics::plot(xtime, mat.na[j, ], xlim = xtimelim, ylim = ylim, 
           xlab = "", ylab = "", axes = FALSE, type = type, 
           col = colour[j], lty = as.numeric(lty[j]), bty="n",
           lwd = as.numeric(lwd[j]), pch=as.numeric(pch[j]))
      graphics::par(new = TRUE)
    }
    ovec$data <- mat.na
    ovec$time <- xtime
    ovec$labs <- labs
  }
  else {
    if (!is.null(labs)) {
      outmat <- NULL
      outlabs <- NULL
      for (j in unique(labs)) {
        temp <- labs == j
        vals <- mat.na[temp, ]
        if (is.matrix(vals)) {
          mvals <- apply(vals, 2, mean, na.rm = TRUE)
        }
        else {
          mvals <- vals
        }
        outmat <- rbind(outmat, mvals)
        outlabs <- c(outlabs, j)
      }
    }
    else {
      outmat <- apply(mat.na, 2, mean, na.rm = TRUE)
      outmat <- rbind(outmat)
      outlabs <- 1
    }
    col.code <- match(col.lty$legend$lab, unique(labs))
    colour <- col.lty$legend$col
    lty <- col.lty$legend$lty
    lwd <- col.lty$legend$lwd
    pch <- col.lty$legend$pch
    for (j in 1:nrow(outmat)) {
      graphics::plot(xtime, outmat[j, ], xlim = xtimelim, ylim = ylim, 
           xlab = "", ylab = "", axes = FALSE, type = type, bty="n", 
           col = colour[col.code[j]], lty = as.numeric(lty[col.code[j]]), 
           lwd = as.numeric(lwd[col.code[j]]), pch = as.numeric(pch[col.code[j]]))
      graphics::par(new = TRUE)
    }
    ovec$data <- outmat
    ovec$time <- xtime
    ovec$labs <- outlabs
  }
  
  if (is.logical(legend)) {
    if (legend) 
    {
      legend <- "topright"
      if ((type=="l") | is.null(pch))
        legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
               lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd))
      else if(type=="p")
        legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
               pch = as.numeric(col.lty$legend$pch) )
      else
        legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
               lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd),pch = as.numeric(col.lty$legend$pch) )
    }
  }
  else 
  {
    if ((type=="l") | is.null(pch))
      legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
             lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd))
    else if(type=="p")
      legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
             pch = as.numeric(col.lty$legend$pch) )
    else
      legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
             lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd),pch = as.numeric(col.lty$legend$pch) )
    
  }
  
  invisible(ovec)
  
}



##' @export
`dplot.norm` <- function (x, labs = NULL, average = FALSE, xlim = NULL, 
                          ylim = NULL,  lty = FALSE, type = "l", colour = TRUE, 
                          lwd = NULL, pch = NULL, legend = "topright",  n = 20) 
{
  ovec <- NULL
  if (is.null(ylim)) 
    ylim <- range(x$data)
  if (is.null(xlim)) 
    xlim <- c(0, 1)
  if (is.null(labs)) {
    labs <- rep(1, nrow(x$index))
  }
  col.lty <- mu.colour(labs, colour, lty, lwd, pch)
  colour <- col.lty$colour
  lty <- col.lty$linetype
  lwd <- col.lty$lwd
  pch <- col.lty$pch
  mat.na <- linear(x, n)
  mat.na$ftime <- x$ftime
  class(mat.na) <- "trackdata"
  xvec <- seq(0, 1, length = n)
  lval <- nrow(x$index)
  if (!average) {
    for (j in 1:lval) {
      graphics::plot(xvec, mat.na[j]$data, xlim = xlim, ylim = ylim, 
           xlab = "", ylab = "", axes = FALSE, type = type, bty="n", 
           col = colour[j], lty = as.numeric(lty[j]), 
           lwd = as.numeric(lwd[j]), pch = as.numeric(pch[j]))
      graphics::par(new = TRUE)
    }
    ovec$data <- mat.na
    ovec$time <- xvec
    ovec$labs <- labs
  }
  else {
    if (!is.null(labs)) {
      outmat <- NULL
      outlabs <- NULL
      for (j in unique(labs)) {
        temp <- labs == j
        vals <- mat.na[temp]$data
        vals <- matrix(vals, ncol = n, byrow = TRUE)
        if (is.matrix(vals)) {
          mvals <- apply(vals, 2, mean)
        }
        else {
          mvals <- vals
        }
        outmat <- rbind(outmat, mvals)
        outlabs <- c(outlabs, j)
      }
    }
    else {
      outmat <- apply(matrix(mat.na, ncol = 20, byrow = TRUE), 
                      2, mean)
      outmat <- rbind(outmat)
      outlabs <- 1
    }
    col.code <- match(col.lty$legend$lab, unique(labs))
    colour <- col.lty$legend$col
    lty <- col.lty$legend$lty
    lwd <- col.lty$legend$lwd
    pch <- col.lty$legend$pch
    for (j in 1:nrow(outmat)) {
      graphics::plot(xvec, outmat[j, ], xlim = xlim, ylim = ylim, 
           xlab = "", ylab = "", axes = FALSE, type = type, bty="n", 
           col = colour[col.code[j]], lty = as.numeric(lty[col.code[j]]), 
           lwd = as.numeric(lwd[col.code[j]]), pch = as.numeric(pch[col.code[j]]))
      graphics::par(new = TRUE)
    }
    ovec$data <- outmat
    ovec$time <- xvec
    ovec$labs <- labs
  }
  
  if (is.logical(legend)) {
    if (legend) 
    {
      legend <- "topright"
      if ((type=="l") | is.null(pch))
        legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
               lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd))
      else if(type=="p")
        legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
               pch = as.numeric(col.lty$legend$pch) )
      else
        legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
               lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd),pch = as.numeric(col.lty$legend$pch) )
    }
  }
  else 
  {
    if ((type=="l") | is.null(pch))
      legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
             lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd))
    else if(type=="p")
      legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
             pch = as.numeric(col.lty$legend$pch) )
    else
      legend(legend, NULL, col.lty$legend$lab, col = col.lty$legend$col, 
             lty = as.numeric(col.lty$legend$lty), lwd = as.numeric(col.lty$legend$lwd),pch = as.numeric(col.lty$legend$pch) )
    
  }
  invisible(ovec)
}

