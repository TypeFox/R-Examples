#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#  GNU General Public License for more details.
#
#  A copy of the GNU General Public License is available at
#  ../../COPYING


###############################################################################
# FUNCTION:                 DESCRIPTION:
#  .xtplot.timeSeries        Plots a 'timeSeries' object
#                            ... support for at = c("pretty", "chic")
###############################################################################
# FUNCTION:                 DESCRIPTION:
# .xtsPlot                   Internal xts plot unitility
# .axTicksByTime2            Takes care of "chic" axis  creation
# .endpoints2                ... determines appropriate axis end points
# .periodicity2              ... determines appropriate axis periodicity
# .colorwheelPalette
###############################################################################


# .plot.timeSeries <-
#   function(
#     x, y, FinCenter = NULL,
#     plot.type = c("multiple", "single"),
#     format = "auto", at = pretty(x),
#     widths = 1, heights = 1,
#     xy.labels, xy.lines, panel = lines, nc, yax.flip = FALSE,
#     mar.multi = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1),
#     oma.multi = c(6, 0, 5, 0), axes = TRUE, ...)
      
#  x=dummySeries();  y = NULL; FinCenter = NULL;  plot.type = "s"    
#  format = "auto"; at = "pretty"; panel = lines; yax.flip = FALSE
#  mar.multi = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1)
#  oma.multi = c(7.75, 1.1, 6.1, 1.1)
#  dots <- list()


.xtplot.timeSeries <-
function(
  x,  y = NULL, FinCenter = NULL,           
  plot.type = c("single", "multiple"),
  format = "auto", at = c("pretty", "chic"),
  panel = lines, yax.flip = FALSE, 
  mar.multi = c(0, 5.1, 0, if (yax.flip) 5.1 else 2.1), 
  oma.multi = c(7.75, 1.1, 6.1, 1.1), # oma.multi = c(6, 0, 5, 0),
  axes = TRUE,
  ...)
{
  # A function implemented by Diethelm Wuertz

  # Description:
  #   Plots timeSeries objects - Internal Function

  # Details:
  #   A modified copy of R's internal 'plotts()' function,
  #   see 'plot.ts()'.

  # FUNCTION:
  
  dots <- list(...)
  if (is.null(dots$minor.ticks)) minor.ticks <- "auto" else minor.ticks <- dots$minor.ticks
  if (is.null(dots$type)) type <- "l" else type <- dots$type
  if (is.null(dots$col)) col <- 1 else col <- dots$col
  if (is.null(dots$pch)) pch <- 20 else pch <- dots$pch
  if (is.null(dots$cex)) cex <- 1 else cex <- dots$cex
  if (is.null(dots$lty)) lty <- 1 else lty <- dots$lty
  if (is.null(dots$lwd)) lwd <- 1 else lwd <- dots$lwd
  if (is.null(dots$grid)) grid <- TRUE else grid <- dots$grid
  if (is.null(dots$frame.plot)) frame.plot <- TRUE else frame.plot <- dots$frame.plot
  if (is.null(dots$ann)) ann <- TRUE else ann <- dots$ann
  if (is.null(dots$cex.axis)) cex.axis <- 1 else cex.axis <- dots$cex.axis
  if (is.null(dots$cex.lab)) cex.lab <- 1 else cex.lab <- dots$cex.lab
  if (is.null(dots$cex.pch)) cex.pch <- 1 else cex.pch <- dots$cex.pch
  
  # Continue ...
  if (minor.ticks == "auto") minor.ticks <- .periodicity2(x)$units
  if (at[1] == "chic") minor.ticks = TRUE
  if (format != "auto") minor.ticks = TRUE
  
  # FinCenter - take care of it:
  if (!is.null(FinCenter)) {
      finCenter(x) <- FinCenter
      if (!missing(y)) finCenter(y) <- FinCenter
      if (is(at, "timeDate")) at@FinCenter <- FinCenter }
  
  # Plot Type:
  plot.type <- plot.type[1]
  if(isUnivariate(x)) plot.type <- "single"
  if(is.timeSeries(y)) plot.type <- "scatter"
  
  # Axis Positions and Format:
  AT <- at[1]
  FORMAT <- format[1]
  if (x@format == "counts") FORMAT <- "counts"
    
  # Main and Labels (We dont use them):
  main <- xlab <- ylab <- ""
  nm <- colnames(x)
  if(length(nm) > 1 &&  ( plot.type == "single" || plot.type == "s"))
    nm <- "Values"

  # Decorations:
  # if (is.null(col)) col <- 1:ncol(x)
  # if (col[1] == 0) col = 1 else col <- .colorwheelPalette(ncol(x)) 
  # if (is.null(pch)) pch <- 20
  # if (is.null(cex)) cex <- 1
  # if (is.null(lty)) lty <- 1
  # if (is.null(lwd)) lwd <- 2
  
  if(is.null(type[1])) type <- "l"
  if (length(type) == 1) type = rep(type, times=NCOL(x))
  if (length(col) == 1) col = rep(col, times=NCOL(x))
  if (length(pch) == 1) pch = rep(pch, times=NCOL(x))
  if (length(cex) == 1) cex = rep(cex, times=NCOL(x))
  if (length(lty) == 1) lty = rep(lty, times=NCOL(x))
  if (length(lwd) == 1) lwd = rep(lwd, times=NCOL(x))
  if (length(cex.pch) == 1) cex.pch = rep(cex.pch, times=NCOL(x))

  X <- as.POSIXct(time(x))
  Y <- series(x)
  
  if (AT == "pretty") {
    at <- pretty(x) }
  if (AT == "chic" ) {
    ep <- .axTicksByTime2(x, format=FORMAT)
    at <- time(x)[ep] } 
  TIME <- time(x)
  
  
  # SINGLE PLOT:
   
  if (plot.type == "single" || plot.type == "s") {

    # All curves in one Frame:
    if (is.null(dots$ylim)) ylim <- range(Y, na.rm = TRUE) else ylim <- dots$ylim
    
    plot(X, Y[, 1], type= "n", ylim = ylim,
     axes = FALSE, main = "", xlab = "", ylab = "")
    for (i in 1:ncol(x)) {
      lines(X, series(x)[, i], type = type[i],
        col = col[i], lty = lty[i], lwd = lwd[i], pch = pch[i],
        cex = cex.pch[i]) 
    }
  
    if (ann) {
      title(main = main, xlab = xlab, ylab = nm, cex.lab = cex.lab) }
  
    if (axes) {
      # X - Axis:
      if (AT == "counts") {
        axis(1, cex.axis = cex.axis)
      } else if (AT == "pretty") {
        
        at <- pretty(time(x))
        if (FORMAT == "auto") format <- "%Y-%m-%d"
        if (!is.null(minor.ticks)) {
          minor.at <- timeSequence(time(x)[1], time(x)[nrow(x)], 
            by = minor.ticks)
          axis.POSIXct(1, at=minor.at, labels=FALSE, col='#BBBBBB', 
            cex.axis = cex.axis) }
        axis.POSIXct(1, at = at, format = format, cex.axis = cex.axis)
        
      } else if (AT == "chic" ) {
        
        ep <- .axTicksByTime2(x, format=FORMAT)
        if (minor.ticks) axis.POSIXct(1, at=TIME, labels=FALSE, col='#BBBBBB',
          cex.axis = cex.axis)
        axis.POSIXct(1, at = TIME[ep], labels=names(ep), 
          las=1, lwd=1, mgp=c(3, 2, 0), cex.axis = cex.axis)
      } else {
        if (minor.ticks) 
              axis.POSIXct(1, at=TIME, labels=FALSE, col='#BBBBBB', 
                cex.axis = cex.axis)
        axis.POSIXct(1, at = at, format = format, cex.axis=cex.axis)
        
      }
      
      # Y - Axis:
      axis(2, cex.axis = cex.axis)
    }
  
    if (frame.plot) {
      box("plot")
    }
  
    if(grid) {
      abline(v = at, lty = 3, col = "darkgrey")
      grid(NA, NULL, lty = 3, col = "darkgrey")
    }
  
    return(invisible())
  }
  
  
  # MULTIPLE PLOT:
  
  if (plot.type == "multiple" || plot.type == "m") {
     
    nser <- ncol(x)
    nc <- if (nser > 4) 2 else 1
    nr <- ceiling(nser/nc)
  
    oldpar <- par(mar = mar.multi, oma = oma.multi, mfcol = c(nr, nc))
    on.exit(par(oldpar))

    for (i in 1:nser) {
        
      plot(X, Y[, i], axes = FALSE, ann = TRUE, type = "n",
        xlab = "", ylab = "",  # log = log, 
        col = col[i], pch = pch[i], lty = lty[i], lwd = lwd[i], cex = cex[i])
      panel(X, Y[, i], type = type[i],
        xlab = "", ylab = "", col = col[i], pch = pch[i], lty = lty[i], 
        lwd = lwd[i], cex = cex.pch[i])
    
      y.side <- if (i%%2 || !yax.flip) 2 else 4
      do.xax <- i%%nr == 0 || i == nser
    
      if (frame.plot) {
        box() }
    
      if (axes) {
        axis(y.side, xpd = NA, cex.axis=cex.axis)
        if (do.xax) {
          if (AT == "counts") {
            axis(1, cex.axis = 1.2 * cex.axis)
          } else if (AT == "pretty") {
            at <- pretty(time(x))
            if (FORMAT == "auto") format <- "%Y-%m-%d"
            TIME <- time(x)
            if (!is.null(minor.ticks)) {
              minor.at <- timeSequence(
                time(x)[1], time(x)[nrow(x)], by=minor.ticks)
              axis.POSIXct(1, at=minor.at, labels=FALSE, 
                cex.axis = 1.2 * cex.axis, col='#BBBBBB') }
            axis.POSIXct(1, at = at, format = format, 
              cex.axis = 1.2 * cex.axis)     
          } else if (AT == "chic" ) {
            ep <- .axTicksByTime2(x, format=FORMAT)
            at <- time(x)[ep]
            format <- attr(ep, "format")
            formatLabels <- names(ep)
            TIME <- time(x)
            if (minor.ticks) 
              axis.POSIXct(1, at=TIME, labels=FALSE, col='#BBBBBB',
                cex.axis = 1.2 * cex.axis)
            axis.POSIXct(1, at = TIME[ep], labels=names(ep), 
              las=1, lwd=1, mgp=c(3, 2, 0), cex.axis = 1.2 * cex.axis)
          } else {
            TIME <- time(x)
            if (minor.ticks) 
              axis.POSIXct(1, at=TIME, labels=FALSE, col='#BBBBBB', 
                cex.axis = 1.2 * cex.axis)
            axis.POSIXct(1, at = at, format = format, 
              cex.axis = 1.2 *cex.axis)
          }
        }
      }
    
      if (ann) {
        mtext(text = nm[i], side = y.side, line = 3, cex = cex.lab)
        if (do.xax) mtext(xlab, side = 1, line = 3, cex = cex.lab)
      }
    
      if(grid) {
        abline(v = at, lty = 3, col = "darkgrey")
        grid(NA, NULL, lty = 3, col = "darkgrey")
      }
          
    } # end of nser loop
  
    return(invisible())
  }
  plot
  
  # SCATTER PLOT:
  
  if (!is.null(y)) {
    
    stopifnot (isUnivariate(x))
    stopifnot (isUnivariate(y))
    plot(series(x), series(y), xlab="", ylab="", col=col, pch=pch, cex=cex)
    
    return(invisible())
    
  }
  
}


###############################################################################
# Test function for xts-plot-like axis positions and labels.


.xtsPlot <- 
  function(x, y=NULL,
    type = "l",
    ann = TRUE,
    axes = TRUE,
    major.ticks = 'auto',
    minor.ticks = TRUE,
    major.format = TRUE,
    grid = TRUE,
    box = TRUE,
    ...) 
{
    # A function written by Diethelm Wuertz
    
    # Descroption:
    #     A simple example to test the xts functions to generate
    #     nice axis positions and Lebels
      
    # Example:
    #     x = 100 * cumulated(LPP2005REC[, 2]); xtsPlot(x)
    
    # Settings:
    
    # time.scale <- periodicity2(x)$scale
    ep <- .axTicksByTime2(x, major.ticks, format.labels=major.format)
    
    # PLOT COORDS:
    xycoords <- xy.coords(time(x), x[, 1])
    
    # RAW PLOT:
    plot(xycoords$x, xycoords$y, type=type, axes=FALSE, ann=FALSE, ...)

    # ADD GRID:
    if (grid) {
      abline(v=xycoords$x[ep], col='grey', lty=3)
      grid(NA, NULL) 
    }

    # ADD AXIS:
    if(axes) {
      if(minor.ticks)
      axis(1, at=xycoords$x, labels=FALSE, col='#BBBBBB')
      axis(1, at=xycoords$x[ep], labels=names(ep), las=1, lwd=1, mgp=c(3,2,0)) 
      axis(2)
    }
    
    # ADD BOX:
    box()
}


# -----------------------------------------------------------------------------
# Borrowed from ...
# Package: xts
# Title: eXtensible Time Series
# Version: 0.9-7
# Date: 2013-06-26
# Author: Jeffrey A. Ryan, Joshua M. Ulrich
# Maintainer: Jeffrey A. Ryan <jeff.a.ryan@gmail.com>
# License: GPL (>= 2)
# URL: http://r-forge.r-project.org/projects/xts/
# Packaged: 2014-01-02 18:00:13 UTC; ripley
# NeedsCompilation: yes
# Repository: CRAN
# Date/Publication: 2014-01-02 19:18:28


.axTicksByTime2 <- 
  function(
    x, ticks.on='auto', k=1, 
    labels=TRUE, format.labels=TRUE, ends=TRUE, gt = 2, lt = 30,
    format = "auto") 
{  
    # A modified function borrowed from the xts-package
    
    # Arguments:
    #   x - a 'timeSerie' Object
    
    # Example:
    #     x = 100 * cumulated(LPP2005REC[, 2]); .axTicksByTime2(x)
    
    tick.opts <- c(
      "years", "months", "weeks", "days", "hours", "minutes", "seconds")
    tick.k.opts <- c(
      10, 5, 2, 1, 6, 1, 1, 1, 4, 2, 1, 30, 15, 1, 1)
    
    if (ticks.on %in% tick.opts) {
      cl <- ticks.on[1]
      ck <- k
    } else {
      tick.opts <- paste(rep(tick.opts, c(4, 2, 1, 1, 3, 3, 1)), tick.k.opts)
      is <- structure(rep(0,length(tick.opts)), .Names = tick.opts)
      
      for(i in 1:length(tick.opts)) {
        y <- strsplit(tick.opts[i], ' ')[[1]]
        ep <- .endpoints2(x, y[1], as.numeric(y[2]))
        is[i] <- length(ep) -1
        if(is[i] > lt) 
        break
      }   
      
      nms <- rev(names(is)[which(is > gt & is < lt)])[1]
      cl <- strsplit(nms, " ")[[1]][1]
      ck <- as.numeric(strsplit(nms, " ")[[1]][2])
    }     
    
    if (is.null(cl)) 
      ep <- NULL else  ep <- .endpoints2(x, cl, ck) 

    if(ends) 
      ep <- ep + c(rep(1,length(ep)-1),0)
    
    if (labels) {
      if(is.logical(format.labels) || is.character(format.labels)) {
        # format by level of time detail, and platform 
        unix <- ifelse(.Platform$OS.type=="unix", TRUE, FALSE)
        time.scale <- .periodicity2(x)$scale
        fmt <- ifelse(unix, '%n%b%n%Y', '%b %Y')
        if (time.scale == "weekly" | time.scale == "daily") 
          fmt <- ifelse(unix, '%b %d%n%Y', '%b %d %Y')
        if (time.scale == "minute" | time.scale == "hourly") 
          fmt <- ifelse(unix, '%b %d%n%H:%M', '%b %d %H:%M')
        if (time.scale == "seconds")
          fmt <- ifelse(unix, '%b %d%n%H:%M:%S', '%b %d %H:%M:%S')
        if(is.character(format.labels)) fmt <- format.labels   
        if (format != "auto") fmt <- format
        names(ep) <- format(time(x)[ep], fmt)
      } else {
        names(ep) <- as.character(time(x)[ep])
      }
    }
    attr(ep, "format") <- fmt

    # Return Value:
    ep  
}


################################################################################


.endpoints2 <-
function (x, on = c("months", "years", "quarters", "weeks", "days", 
    "hours", "minutes", "seconds"), k = 1) 
{
    # A modified function borrowed from the xts-package
  
    # Arguments:
    #    x - a 'timeDate' object
  
    # Example:
    #     x <- 100 * cumulated(LPP2005REC[, 2]); .endpoints2(x)
    
    stopifnot(is(x, "timeSeries"))
    x <- time(x)
  
    on <- match.arg(on)
    posix <- as.POSIXct(x, origin="1970-01-01")
    .posix <- unclass(posix)
    if (on == "years") {
        ans <- as.integer(which(diff(as.POSIXlt(posix)$year%/%k + 1) != 0))
    }
    else if (on == "quarters") {
        ans <- as.integer(which(diff((as.POSIXlt(posix)$mon%/%3) + 1) != 0))
    }
    else if (on == "months") {
        ans <- as.integer(which(diff(as.POSIXlt(posix)$mon%/%k + 1) != 0))
    }
    else if (on == "weeks") {
        ans <- as.integer(
          which(diff((.posix + (3L * 86400L))%/%604800L%/%k + 1) != 0))
    }
    else if (on == "days") {
        ans <- as.integer(which(diff(.posix%/%86400L%/%k + 1) != 0))
    }
    else if (on == "hours") {
        ans <- as.integer(which(diff(.posix%/%3600L%/%k + 1) != 0))
    }
    else if (on == "minutes" || on == "mins") {
        ans <- as.integer(which(diff(.posix%/%60L%/%k + 1) != 0))
    }
    else if (on == "seconds" || on == "secs") {
        ans <- as.integer(which(diff(.posix%/%k + 1) != 0))
    }
    ans <- c(0, ans, NROW(x))
  
    # Return Value:
    ans
}


###############################################################################


.periodicity2 <-
   function (x) 
{
    # A modified function borrowed from the xts-package
     
    # Arguments:
    #    x - a 'timeDate' object
  
    # Example:
    #     x <- 100 * cumulated(LPP2005REC[, 2]); .periodicity2(x)
     
    # FUNCTION:
    # Check Argument:
    stopifnot(is(x, "timeSeries"))
    x <- time(x)
    
    p <- median(diff(as.integer(as.POSIXct(x, origin="1970-01-01"))))
     
    if (is.na(p)) stop("cannot calculate periodicity of 1 observation")
    units <- "days"
    scale <- "yearly"
    label <- "year"
     
    if (p < 60) {
        units <- "secs"
        scale <- "seconds"
        label <- "second"
    } else if (p < 3600) {
        units <- "mins"
        scale <- "minute"
        label <- "minute"
        p <- p/60L
    } else if (p < 86400) {
        units <- "hours"
        scale <- "hourly"
        label <- "hour"
    } else if (p == 86400) {
        scale <- "daily"
        label <- "day"
    } else if (p <= 604800) {
        scale <- "weekly"
        label <- "week"
    } else if (p <= 2678400) {
        scale <- "monthly"
        label <- "month"
    } else if (p <= 7948800) {
        scale <- "quarterly"
        label <- "quarter"
    }
     
    # Return Value:
    list(
      difftime = structure(p, units = units, class = "difftime"), 
      frequency = p, start = start(x), end = end(x), units = units, 
      scale = scale, label = label)
}


###############################################################################


.colorwheelPalette <- 
  function(n)
{
  # A function implemented by Diethelm Wuertz
    
  # FUNCTION:
  
  # Color Wheel:
  orig <- c(
    "#FFF200", "#FBAA19", "#F26522", "#EF4823",
    "#ED1D24", "#A9285F", "#662D91", "#4D2F91",
    "#2E3092", "#00707E", "#00A650", "#8CC63F")
  orig <- orig[-1]

  # Splice Wheel
  if (n == 11) return(orig)
    rgb.tim <- t(col2rgb(orig))
    temp <- matrix(NA, ncol = 3, nrow = n)
    x <- seq(0, 1, , 11)
    xg <- seq(0, 1, , n)
    for (k in 1:3) {
        hold = spline(x, rgb.tim[, k], n = n)$y
        hold[hold < 0] = 0
        hold[hold > 255] = 255
        temp[, k] = round(hold)
    }
    ans <- rgb(temp[, 1], temp[, 2], temp[, 3], maxColorValue = 255)
    
    # Return Value:
    ans
}


###############################################################################
 
