tierChart <- function(x,
                      startMonth = latestJanuary(end(x)),
                      nMonths = 4,
                      nYears = 7,
                      offsets = 0,
                      padDays = 6,
                      pch = "year", 
                      lty = "solid",
                      lwd = 1.5, 
                      col = 1 + (n.y:1),
                      type = "b", 
                      ylim = NULL,
                      outlier.trim = 0,
                      noTrimLastYear = TRUE,
                      extendHorizontalTicks = TRUE,
                      circles.ymd = NULL,
                      circles.col = 6,
                      circles.inches = 0.1,
                      vlines.ymd = NULL,
                      vlines.col = 2,
                      vlines.lty = 4,
                      vlines.lwd = 1.5,
                      vlines2.ymd = NULL,
                      vlines2.col = 3,
                      vlines2.lty = "solid",
                      vlines2.lwd = 2,
                      hlines = NULL,
                      hlines.col = 1,
                      hlines.lty = 1,
                      hlines.lwd = 1,
                      tiPoints.1 = NULL,
                      tiPoints.2 = NULL,
                      pch.1 = "*", 
                      pch.2 = "+",
                      col.1 = 2, 
                      col.2 = 3,
                      nolegend = FALSE,
                      main = deparse(substitute(x)),
                      topleft.labels = NULL,
                      topright.labels = NULL,
                      legend.ncol = length(years),
                      legend.bg = 0,
                      timestamp = TRUE,
                      topline = TRUE,
                      vlines.periodEnd = TRUE,
                      vlines.month  = TRUE,
                      midperiod = FALSE,
                      lwdLastYear = 1.5,
                      cex = 1.5,
                      boxes = TRUE,
                      ...){
  ## vlines.monday has been replaced by vlines.periodEnd, but it remains a
  ## valid parameter name for the sake of backwards compatibility
  
  ## graphics parameters
  oldpar <-  par(xpd = FALSE, err = -1, mar = c(2,5,3,5)+0.1, mgp = c(1,0,0),
                 xaxs = "i", yaxs = "i", mex = 1, ann = FALSE)
  on.exit(par(oldpar), add = TRUE) 
  if(par("new") == FALSE) plot.new()
  main <- main
  ## The data lines and points for tier charts are drawn by the matplot()
  ## function. We start by creating  matrices with columns for each
  ## year plotted. We'll need:
  ##   py -- y coordinates of the points to be plotted
  ##   px -- x coordinates of the points to be plotted
  ##   pymd -- ymd dates corresponding to the plotted points
  ##   pindex -- subscripts for x, i.e., px[,i] == x[pindex[,i]]  
  ##             (except for NA values in pindex)
  ## Then we'll draw the plot and return.
  
  x  <- as.tis(x)
  xTif <- tif(x)
  dim(x) <- NULL
  if(!is.ti(startMonth)) stop("startMonth must be a ti")
  xEndMonth <- currentMonth(end(x))
  if(startMonth > xEndMonth) stop("startMonth after end(x)")
  endMonth <- startMonth + nMonths - 1
  if(endMonth > xEndMonth){
    xExtendedEnd <- ti(endMonth, tif(x))
    x[xExtendedEnd] <- NA
  }
  xStartMonth <- currentMonth(start(x))
  
  years  <- (year(startMonth) - nYears + 1): year(startMonth)
  eyears <- (year(endMonth)   - nYears + 1): year(endMonth)
  
  eym <- 100*eyears + month(endMonth)
  years <- years[eym >= (ymd(start(x)) %/% 100)]

  ## julian dates of the series
  if(midperiod) xj <- floor(jul(ti(x), offset = 0.5))
  else          xj <- jul(ti(x))

  ## some lengths and sequence vectors
  n.y      <- length(years)
  n.x      <- length(x)
  y.seq    <- 1:n.y
  x.seq    <- 1:n.x

  startMonths <- startMonth - 12*((nYears - 1):0)
  endMonths   <- startMonths + nMonths - 1
  ## drop (start,end) pairs outside the range of data
  inRange <- startMonths <= xEndMonth & endMonths >= xStartMonth
  startMonths <- startMonths[inRange]
  endMonths   <- endMonths[inRange]
  ## starts and ends: vectors of starting and ending Julian dates for each year
  ## on the plot
  offsets <- c(rev(offsets), rep(0, length(years) - length(offsets)))
  starts <- jul(startMonths - 1) - (padDays - 1) + offsets
  ends   <- jul(endMonths) + padDays + offsets

  ## the x-axis starts at 1 and has daily ticks
  x.ticks <- 1:(1 + ends[n.y] - starts[n.y])
  ## year[n.y] ymd dates associated with x.ticks
  x.ymd <- ymd(seq(starts[n.y], ends[n.y]))
  ## x.periodEnds[i] is TRUE iff x.ymd[i] is the last day of a period of x
  x.periodEnds <- x.ymd == ymd(ti(x.ymd, tif = xTif))
  
  ## If Feb 28 appears on the plot, Feb 29 has to as well. If it's not already
  ## there, add a fake observation to x.ticks, x.ymd and x.periodEnds to
  ## accomodate it. 
  feb28 <- match(0228, x.ymd %% 10000, nomatch = FALSE)
  if(feb28 != FALSE){
    feb29 <- match(0229, x.ymd %% 10000, nomatch = FALSE)
    if(feb29 == FALSE){
      x.ticks       <- append(x.ticks, length(x.ticks) + 1)
      x.ymd         <- append(x.ymd, x.ymd[feb28]+1, after = feb28)
      x.periodEnds  <- append(x.periodEnds, FALSE, after = feb28)
    }
  }
  xtlen <- length(x.ticks)
  
  ## initialize matrices with NA's 
  matrows <- ceiling(xtlen*frequency(x)/365) + 1
  py <- matrix(NA, nrow = matrows, ncol =  n.y,
               dimnames = list(character(0), paste(years)))  
  pindex <- pymd  <- px <- py
  
  ## fill in matrices by looping thru the years
  for(y in y.seq){ 
    ix <- x.seq[between(xj, starts[y], ends[y])]
    iRows <- 1:length(ix)  ## Rows to fill in

    pindex[iRows, y] <- ix
        py[iRows, y] <- x[ix]
      pymd[iRows, y] <- ymd(xj[ix])
          ymd.offset <- 10000 * (max(years) - years[y])
        px[iRows, y] <- match(ymd(xj[ix] + offsets[y]) + ymd.offset, x.ymd)
  }
  ## The plot is drawn in six steps. The first is to draw a blank plot, solely
  ## for the purpose of setting the plot coordinates.  Next, we add the legends
  ## and titles, and only then do we actually draw the data lines and
  ## symbols. Reversing the usual order this way allows the data lines to
  ## overwrite the legends and titles. The fourth step is to draw the grid lines
  ## and the axes. Fifth is adding in the optional "extra" points and vertical
  ## lines specified by tiPoints.1, tiPoints.2, and vlines.ymd. The sixth and
  ## final step is to draw circles at the points indicated by circles.ymd

  ## Step 1: Setting coordinates
  xlim <- c(x.ticks[1], x.ticks[xtlen] + 1)
  ##   The y coordinates are governed by two parameters: outlier.trim and ylim.
  ##   If ylims are set explicitly, outlier.trim has no effect.  Otherwise, ylim
  ##   will be calculated by making sure that the 
  ##   c(outlier.trim, 1 - outlier.trim) quantiles fall within ylim.  The
  ##   default value (outlier.trim = 0) makes sure all observations fall within
  ##   the plot area.
  ##   In either case, ylim will be further stretched (if necessary) to
  ##   accomodate "pretty" y-axis labels.
  if(missing(ylim)){
    bounds <- 1.1 * quantile(py, c(outlier.trim, 1-outlier.trim), na.rm=TRUE) - 
              0.1 * median(py, na.rm=TRUE)
    ylim <- c(max(bounds[1], min(py, na.rm=TRUE)), 
              min(bounds[2], max(py, na.rm=TRUE)))
    if(noTrimLastYear){
      pyLastYear <- py[,ncol(py)]
      ylim <- c(min(bounds[1], min(pyLastYear, na.rm = TRUE)),
                max(bounds[2], max(pyLastYear, na.rm = TRUE)))
    }
  }
  ylim[2] <- 1.1*ylim[2] - 0.1*ylim[1]
  y.ticks <- pretty(ylim, n = 8)
  ylim <- range(y.ticks)
  ##   draw the blank plot.
  plot.window(xlim = xlim, ylim = ylim)

  
  ## Step 2: Add legends and titles
  if( length(pch) == 1 && is.character(pch) && pch == "year")
    pch <- paste(substring(years, 4, 4), collapse = "")

  if(!missing(lwdLastYear)){
    lwd <- rep(lwd, length.out = length(years))
    lwd[length(lwd)] <- lwdLastYear
  }
  
  if(!nolegend){
    legendStrings <- format(years)
    if(any(offsets != 0)){
      roffStrings <- paste("[", offsets, "]", sep = "")
      roffStrings[offsets == 0] <- ""
      legendStrings <- paste(legendStrings, roffStrings, sep = "")
    }
    legend(x = 0.5*xlim[1] + 0.5*xlim[2],
           y = max(ylim),
           legend = legendStrings, 
           pch = pch,
           lty = lty,
           lwd = lwd,
           col = col,
           ncol = legend.ncol,
           bg = legend.bg, bty = "o", cex = 1.3,
           xjust = 0.5)  
  }

  if(length(main)) 
	mtext(main, side = 3, adj = 0.5, cex = 1.5, line = 0.5)	
  ##   upper left labels
  if(tl.len <- length(topleft.labels)){
    for(i in 1:tl.len)
      mtext(topleft.labels[i], side = 3, adj = 0, at = xlim[1], 
	    line = tl.len - i + 0.2)
  }
  ##   upper right labels
  if(tr.len <- length(topright.labels)){
    for(i in 1:tr.len)
      mtext(topright.labels[i], side = 3, adj = 1, at = xlim[2], 
	    line = tr.len - i + 0.2)
  }
  
  ## Step 3: Data lines and points
  par(new = TRUE) ## don't erase what we've already drawn
  matplot(px, py, type = type, pch = pch, col = col, cex = cex, 
          lty = lty, lwd = lwd, axes = FALSE, bty = "n", yaxs = "i", xaxs = "i", 
          xlim = xlim, ylim = range(y.ticks), ...) 
  if(type == "b" | type == "p") 
    matpoints(px, py, pch = pch, lty = 1, col = col, cex = cex)
  
  ## Step 4: Grid lines and axes
  ##   locations of period ends, and days 1 and 15 of a month
  periodEnd.ticks  <- x.ticks[x.periodEnds]
  day15 <- (1:xtlen)[x.ymd %% 100 == 15]
  eom   <- (1:xtlen)[x.ymd == ymd(ti(x.ymd, tif = "monthly"))]
  ##   vertical periodEnd dashed lines
  if(vlines.periodEnd)
    axis(1, at = periodEnd.ticks, tck = 1, lty = 2, labels = FALSE)
  ##   vertical last day of month grid lines
  if(vlines.month) abline(v = eom, lty = 1)
  ##   left axis
  if(extendHorizontalTicks)
    axis(2, at = y.ticks, tck = 1, lty = 2, labels = FALSE) 
  axis(2, at = y.ticks, las = 1, adj = 1, tck = 0.02, labels = TRUE) 
  ##   right axis
  axis(4, at = y.ticks,  las = 1, adj = 0, tck = 0.02) 
  ##   top axis
  if(topline) abline(h = ylim[2])
  ##   bottom axis
  axis(1, at = x.ticks, tck = 0.02, labels = FALSE)

  ##   if the plot covers 8 or more months, use abbreviated month labels
  ##   only. Otherwise use full month labels and PeriodEnd labels
  monthLabels <- format(jul(x.ymd[day15]), "%b %Y")
  if(length(monthLabels) > 8)
    axis(1, at = periodEnd.ticks,  tck = 0.04, labels = FALSE)
  else {
    periodEndLabels <- paste(x.ymd[x.periodEnds] %% 100)
    axis(1, at = periodEnd.ticks,  tck = 0.04, labels = periodEndLabels)
  }
  mtext(text = monthLabels, side = 1, line = 1, at = day15) 

  ## Step 5: "Extra" points and vlines
  if(!is.null(vlines.ymd)){  
    ## draw vertical lines at the dates in vlines.ymd
    i <- match(vlines.ymd, x.ymd, nomatch = 0)
    if(any(i != 0))
      abline(v  =  x.ticks[i], col = vlines.col, lty = vlines.lty, lwd = vlines.lwd)
  }
  if(!is.null(vlines2.ymd)){  
    ## draw vertical lines at the dates in vlines.ymd
    i <- match(vlines2.ymd, x.ymd, nomatch = 0)
    if(any(i != 0))
      abline(v  =  x.ticks[i], col = vlines2.col, lty = vlines2.lty, lwd = vlines2.lwd)
  }
  if(!is.null(hlines)){
    abline(h = hlines, col = hlines.col, lty = hlines.lty, lwd = hlines.lwd)
  }
  if(!is.null(tiPoints.1)){  
    ## add points to plot
    ymd.1 <- ymd(ti(tiPoints.1))
    i.1 <- match(ymd.1, x.ymd, nomatch = 0)
    if(any(i.1 != 0))
      matpoints(x.ticks[i.1], tiPoints.1[i.1 > 0], pch = pch.1,
		cex = 1.4, col = col.1)
  }

  if(!is.null(tiPoints.2)){  
    ## add even more points to plot
    ymd.2 <- ymd(ti(tiPoints.2))
    i.2 <- match(ymd.2, x.ymd, nomatch = 0)
    if(any(i.2 != 0))
      matpoints(x.ticks[i.2], tiPoints.2[i.2 > 0], pch = pch.2,
		cex = 1.4, col = col.2)
  }
  
  ## Step 6 : Circle the points whose dates match circles.ymd
  if(!is.null(circles.ymd)){
    i.circles <- match(circles.ymd, pymd, nomatch = 0)
    if(any(i.circles != 0)){
      x <- px[i.circles]
      y <- py[i.circles]
      symbols(x, y, add = TRUE, circles = rep(1, length(x)),
	      inches = circles.inches, col = circles.col)
    }
  }

  if(boxes){
    ## add arrow boxes 0.5 x 0.25 inches
    extents <- c(diff(xlim), diff(ylim)) / (par("pin") * c(2, 4))
    boxBottom <- ylim[2] - extents[2]
    lBox <- c(xlim[1], boxBottom, xlim[1]+extents[1], ylim[2])
    rBox <- c(xlim[2]-extents[1], boxBottom, xlim[2], ylim[2])
    rect(lBox[1], lBox[2], lBox[3], lBox[4],  col = "white")
    rect(rBox[1], rBox[2], rBox[3], rBox[4],  col = "white")
    arrows(lBox[3], (lBox[2] +lBox[4])/2, lBox[1], (lBox[2]+lBox[4])/2,
           angle = 30,lwd = 4, col = "black")
    arrows(rBox[1], (rBox[2] +rBox[4])/2, rBox[3], (rBox[2]+rBox[4])/2,
           angle = 30,lwd = 4, col = "black")
    printBox <- c(lBox[3], lBox[2], lBox[3] + (lBox[3] - lBox[1]), lBox[4])
    rect(printBox[1], printBox[2], printBox[3], printBox[4], col = "darkblue")
    text(x = mean(printBox[c(1,3)]),
         y = mean(printBox[c(2,4)]),
         labels = "Print Me",
         adj = c(0.5, 0.5),
         cex = 1.25, col = "white")
  }
  else {
    lBox <- rBox <- printBox <- NULL
  }
  
  ## Date-Time Stamp 
  if(timestamp) mtext(date(), side = 3, outer = TRUE, line = -1, adj = 1)
  
  invisible(list(px = px, py = py, ymd = pymd, index = pindex,
                 lBox = lBox, rBox = rBox, printBox = printBox,
                 startMonth = startMonth, nMonths = nMonths,
                 nYears = nYears))
}

adjustableTierChart <- function(x, ..., edit = TRUE,
                                changes = numeric(0), verbose = FALSE){
  locateOne <- function(){
    options(error = NULL)
    point <- try(locator(1), silent = TRUE)
    if(class(point) == "try-error")
      point <- NULL
    options(error = quote(recover()))
    point
  }

  arglist <- list(..., edit = edit, verbose = verbose)
  if(is.null(arglist$main)) 
    arglist$main <- deparse(substitute(x))
  arglist$x <- x
  arglist$changes <- changes
  arglist$boxes <- TRUE
  drops <- match(c("edit", "changes", "verbose"), names(arglist), nomatch = 0)
  tcArgs <- arglist[-drops]
  tc <- do.call("tierChart", tcArgs)

  xStart <- start(x)
  xEnd   <- end(x)
  xStartMonth <- ti(xStart, tif = "monthly")
  xEndMonth   <- ti(xEnd,   tif = "monthly")
  
  ## adjust
  lastcol    <- dim(tc$px)[2]
  not.na     <- !(is.na(tc$px[,lastcol]) | is.na(tc$py[,lastcol]))
  px         <- tc$px[not.na, lastcol]
  py         <- tc$py[not.na, lastcol]
  ymd        <- tc$ymd[not.na, lastcol]
  index      <- tc$index[not.na, lastcol]
  xStep      <- mean(diff(px))
  lBox       <- tc$lBox
  rBox       <- tc$rBox
  pBox       <- tc$printBox
  startMonth <- tc$startMonth
  nMonths    <- tc$nMonths
  scrollIncrement <- max(nMonths - 2, 1)
  if(frequency(xStart) <= 12) scrollIncrement <- round(frequency(xStart))
  
  ## refresh arrows for changes that are still in range
  if(length(changes) > 0){
    for(i in 1:nrow(changes)){
      change <- changes[i,]
      j <- match(change["ymd"], ymd, nomatch = 0)
      if(j > 0) arrows(px[j], py[j], px[j], change["newY"], col = 2)
    } 
  }
  
  if(verbose)  cat("Awaiting mouse input...")
  nextpoint  <- locateOne()
    
  while(length(nextpoint) != 0){
    npX <- nextpoint$x
    npY <- nextpoint$y
    lBoxClicked <- rBoxClicked <- pBoxClicked <- FALSE
    ## see if user clicked in lBox, rBox or pBox
    if((npY >= lBox[2]) && (npY <= lBox[4])){
      lBoxClicked <- (npX >= lBox[1]) && (npX <= lBox[3])
      rBoxClicked <- (npX >= rBox[1]) && (npX <= rBox[3])
      pBoxClicked <- (npX >= pBox[1]) && (npX <= pBox[3])
    }
    if(rBoxClicked){
      inc <- min(scrollIncrement, xEndMonth - startMonth)
      arglist$startMonth <- startMonth + inc
      return(do.call("adjustableTierChart", arglist))
    }      
    if(lBoxClicked){
      inc <- min(scrollIncrement, startMonth - xStartMonth)
      arglist$startMonth <- startMonth - inc
      return(do.call("adjustableTierChart", arglist))
    }
    if(pBoxClicked){
      printer <- options()$colorPrinter
      if(is.null(printer)) cat("options()$colorPrinter not set\n")
      else {
        postscript(file = "",
                   command = paste("lp -d", printer, sep = ""),
                   pointsize = 9)
        setColors(frColors())
        modArgs <- tcArgs
        modArgs$legend.bg <- "grey90"
        do.call("tierChart", modArgs)
        dev.off()
      }
    }
    else {
      if(edit){
        ## find nearest point (by x value)
        xDistances <- abs(px - npX)
        if(min(xDistances) < xStep){ 
          i <- order(xDistances)[1]
          ## change the yvalue
          new.y <- npY
          ## draw an arrow from old point to new point,
          ##   color it the same color 
          ##   as the most recent year plotted
          changes <- rbind(changes, c(ymd = ymd[i], newY = npY))
          arglist$changes <- changes
          arrows(px[i], py[i], px[i], new.y, col = 2)
          if(verbose)
            cat("Move", ymd[i], "obs from", signif(py[i], 3),
                "to", signif(new.y, 3), "\nAwaiting mouse input... ")
        }
      }
    }
    ## repeat
    nextpoint <- locateOne()
  }

  if(length(changes) > 0)
    x[match(changes[,"ymd"], ymd(ti(x)))] <- changes[,"newY"]
  retval <- x
  attr(retval, "startMonth") <- arglist$startMonth
  if(verbose) cat("adjustableTierChart() exiting.\n")
  invisible(retval) 
}

