scatterPlot <- function(x, y, 
                        plotType = "p",
                        lineType = "solid",
                        lineWidth = 1.5,
                        plotChar = "*",
                        dataCex = 1,
                        color = "black",
                        ## xAxisSpec
                        xAxisMin = NULL,
                        xAxisMax = NULL,
                        xExpandBy = 0.04,
                        xTicks = 5,
                        xTickLocations = NULL,
                        labelXTicks = TRUE,
                        xTickLabels = NULL,
                        xCex = 1,
                        xAxisLabel = NULL,
                        labelXAxis = TRUE,
                        xSpace = 4,
                        ## yAxisSpec
                        yAxisMin = NULL,
                        yAxisMax = NULL,
                        yExpandBy = 0.04,
                        yTicks = 5,
                        yTickLocations = NULL,
                        yTickLabels = NULL,
                        labelLeftTicks = FALSE,
                        labelRightTicks = TRUE,
                        yCex = 1,
                        extendTopTick = TRUE,
                        leftAxisLabel = NULL,
                        rightAxisLabel = NULL,
                        labelLeftAxis = TRUE,
                        labelRightAxis = FALSE,
                        ## textItemsSpec
                        cex = 1,
                        head = NULL,
                        headAlign = 0.5,
                        headCex = 1.5,
                        sub = NULL,
                        subCex = 0.85,
                        leftTopLabel = NULL,
                        rightTopLabel = NULL,
                        topLabelAlign = 0,
                        labCex = 1,
                        leftInsideLabel = NULL,
                        rightInsideLabel = NULL,
                        innerOffset = 0.05,
                        innerCex = 0.8,
                        foot = NULL,
                        footAlign = 0,
                        footCex = 0.8,
                        footSpace = -1,
                        ## layoutSpec
                        tck = 0.03,
                        axisWidth = 2,
                        boxType = "u",
                        leftMargin = -1,
                        rightMargin = -1,
                        topMargin = -1,
                        bottomMargin = -1){
  if(missing(cex)) cex <- par("cex")
  ## Reset some graphical parameters
  par.old <- par(xaxs = "i", yaxs = "i", mex = 1)
  if(par("new") == FALSE)
    plot.new()

  if(missing(x) && missing(y)){
    ## create dummy arguments and set color to the background color
    x <- y <- c(0, 1)
    xyNames <- c("dummyX", "dummyY")
    color <- 0
  }
  else {
    if(missing(x) || missing(y))
      stop("Must supply x and y args to scatterPlot()")
    xyNames <- c(deparse(substitute(x)), deparse(substitute(y)))
    if(is.character(x)){
      xyNames[1] <- x
      x <- evalOrEcho(x)
    }
    if(is.character(y)){
      xyNames[2] <- y
      y <- evalOrEcho(y)
    }
  }

  ## use intersection of x and y
  if(is.tis(x) && is.tis(y)) xy <- cbind(x,y)
  else                     xy <- cbind(x, y)
  
  colnames(xy) <- xyNames
  x <- xy[,1]
  y <- xy[,2]

  ## Scape specific cex's by overall cex
  dataCex  <- cex * dataCex
  headCex  <- cex * headCex
  subCex   <- cex * subCex
  labCex   <- cex * labCex
  footCex  <- cex * footCex
  xCex     <- cex * xCex
  yCex     <- cex * yCex
  innerCex <- cex * innerCex

  ## Plotting is done in the following order:
  ##  the figure margins (axes, titles, footnotes)
  ##  the lines and/or points representing the data
  ##  legends
  expandRange <- function(inRange, expandBy = 0.04){
    ## local function
    if(length(inRange) != 2) stop("length(inRange) != 2")
    addFactor <- expandBy * (inRange[2] - inRange[1])
    return(inRange + c(-addFactor, addFactor))
  }

  ## set xRange and yRange
  xDataRange <- range(x)
  xRange <- expandRange(xDataRange, xExpandBy)
  if(!is.null(xAxisMin)) xRange[1] <- xAxisMin
  if(!is.null(xAxisMax)) xRange[2] <- xAxisMax
  
  yDataRange <- range(y)
  yRange <- expandRange(yDataRange, yExpandBy)
  if(!is.null(yAxisMin)) yRange[1] <- yAxisMin
  if(!is.null(yAxisMax)) yRange[2] <- yAxisMax
  
  yAt <- NULL
  ## figure out yAt based on yTicks and yTickLocations
  if(!is.null(yTickLocations))
    yAt <- evalOrEcho(yTickLocations)
  if(is.null(yAt)){
    if(yTicks > 0){
      yAt <- signif(pretty(yRange, n = yTicks))
      yAt <- yAt[yAt >= yRange[1]]
    }
  }

  if(extendTopTick && !is.null(yAt))
    yRange[2] <- max(c(yRange[2], yAt))

  ## Set tick labels
  if(labelLeftTicks || labelRightTicks){
    if(is.null(yTickLabels) && !is.null(yAt))
      yTickLabels <- format(signif(yAt, 5))
  }
  else yTickLabels <- FALSE

  xAt <- NULL
  ## figure out xAt based on xTicks and xTickLocations
  if(!is.null(xTickLocations))
    xAt <- evalOrEcho(xTickLocations)
  if(is.null(xAt)){
    if(xTicks > 0){
      xAt <- signif(pretty(xRange, n = xTicks))
      xAt <- xAt[xAt >= xRange[1]]
    }
  }
  if(labelXTicks){
    if(is.null(xTickLabels) && !is.null(xAt))
      xTickLabels <- format(signif(xAt, 5))
  }
  else xTickLabels <- FALSE

  ## By default, label only one of left and right axes
  if(missing(labelRightAxis) && labelLeftAxis)  labelRightAxis <- FALSE
  if(missing(labelLeftAxis)  && labelRightAxis) labelLeftAxis  <- FALSE
  
  ## FIGURE MARGINS -- includes titles, footnotes, and axes
  ## Before these can be plotted, we have to calculate the required
  ## margin sizes and set par("mar") appropriately.

  ## The top figure margin needs space for head, sub, leftTopLabel, and
  ## rightTopLabel.

  ## Before we can do anything with text items, we have to "decode" them.
  ## Text item arguments may be given as simple strings, collections of
  ## strings, or as expressions meant to be evaluated to create strings.  To
  ## decode the item, the function evalOrEcho() is used.  evalOrEcho attempts
  ## to parse and evaluate its argument, and return the result.  If it cannot
  ## do so, it simply returns its argument.  So expressions get evaluated,
  ## while simple strings just get echoed.  All of these forms work:
  ##
  ## 1.   sub = "This is a simple subtitle"
  ## 2.   sub = c("this is a two", "line subtitle")
  ## 3.   sub = 'c("this is another", "two line subtitle")'
  ##
  head             <- evalOrEcho(head, resultMode = "character")
  sub              <- evalOrEcho(sub, resultMode = "character")
  foot             <- evalOrEcho(foot, resultMode = "character")
  leftTopLabel     <- evalOrEcho(leftTopLabel, resultMode = "character")
  rightTopLabel    <- evalOrEcho(rightTopLabel, resultMode = "character")
  leftInsideLabel  <- evalOrEcho(leftInsideLabel, resultMode = "character")
  rightInsideLabel <- evalOrEcho(rightInsideLabel, resultMode = "character")
  xAxisLabel       <- evalOrEcho(xAxisLabel, resultMode = "character")
  leftAxisLabel    <- evalOrEcho(leftAxisLabel, resultMode = "character")
  rightAxisLabel   <- evalOrEcho(rightAxisLabel, resultMode = "character")

  if(labelXAxis     && is.null(xAxisLabel))         xAxisLabel <- xyNames[1]
  if(labelLeftAxis  && is.null(leftAxisLabel))   leftAxisLabel <- xyNames[2]
  if(labelRightAxis && is.null(rightAxisLabel)) rightAxisLabel <- xyNames[2]
  
  nHeads             <- length(head)
  nSubs              <- length(sub)
  nLeftTopLabels     <- length(leftTopLabel)
  nRightTopLabels    <- length(rightTopLabel)
  nLeftInsideLabels  <- length(leftInsideLabel)
  nRightInsideLabels <- length(rightInsideLabel)
  nXAxisLabels       <- length(xAxisLabel)
  nLeftAxisLabels    <- length(leftAxisLabel)
  nRightAxisLabels   <- length(rightAxisLabel)

  ylabSpace <- labCex * max(c(nLeftTopLabels, nRightTopLabels))
  headSpace <- headCex * nHeads
  subSpace <- subCex * nSubs
  topSpace <- max(c(ylabSpace, headSpace + subSpace))
  
  ## The bottom figure margin gets 4 lines of space for ticklabels,
  ## yearLabels, and/or xAxisLabel, plus whatever space is needed for
  ## footnotes.
  if(footSpace < 0) footSpace <- length(foot)
  bottomSpace <- xSpace * labCex + footSpace * footCex

  if(is.character(yTickLabels))
    tickLabelSpace <- max(strwidth(yTickLabels, units = "inches", cex = yCex)/
                          (par("csi")*par("mex")))
  else
    tickLabelSpace <- 0

  leftTickLabelSpace  <- tickLabelSpace * labelLeftTicks
  rightTickLabelSpace <- tickLabelSpace * labelRightTicks
  leftSpace  <- 1 + leftTickLabelSpace  + nLeftAxisLabels
  rightSpace <- 1 + rightTickLabelSpace + nRightAxisLabels
  
  mar <- par("mar")
  if( leftMargin  < 0)   leftMargin <- max( leftSpace, mar[2])
  if(rightMargin  < 0)  rightMargin <- max(rightSpace, mar[4])
  if(   topMargin < 0)    topMargin <- max(  topSpace, 1)
  if(bottomMargin < 0) bottomMargin <- bottomSpace

  par(mar = c(bottomMargin, leftMargin, topMargin, rightMargin))

  ## Now that margins have been set, set the plot window's user coords
  plot.window(xRange, yRange)

  ## Place titles and footnotes
  ##   head
  if(nHeads)
    for(i in 1:nHeads)
      mtext(head[i], side = 3, adj = headAlign, cex = headCex,
            line = headSpace + subSpace - i * headCex)
  ##   sub
  if(nSubs)
    for(i in 1:nSubs)
      mtext(sub[i], side = 3, adj = headAlign, cex = subCex,
            line = subSpace - i * subCex)
  ##   foot
  if(len.foot <- length(foot))
    for(i in 1:len.foot)
      mtext(foot[i], side = 1, adj = footAlign, cex = footCex,
            line = xSpace * labCex + (i - 1) * footCex)
  
  ## plot the points
  plot.xy(xy.coords(x, y), type = plotType, lty = lineType,
          lwd = lineWidth, col = color, pch = plotChar, cex = dataCex) 

  box(lwd = axisWidth, bty = boxType)

  ## draw the X axis.
  
  ## mgp parameter is supposed to be the margin line in mex units for
  ## c(axis title, axis labels, axisline).  We don't use the first, and the
  ## last is always 0, but we need to compute a reasonable number for the
  ## second.
  xmgp <- ymgp <- par("mgp")
  xmgp[2] <- max(0.6*xCex - 0.5, 0)
  ymgp[2] <- max(0.6*yCex - 0.2, 0)
  
  axis(side = 1, at = xAt, labels = if(labelXTicks) xTickLabels else FALSE,
       lwd = axisWidth, cex.axis = xCex, tick = TRUE, tck = tck,
       mgp = xmgp)
  ## draw left and right Axes
  axis(side = 2, at = yAt, labels = if(labelLeftTicks) yTickLabels else FALSE,
       lwd = axisWidth, cex.axis = yCex, tick = TRUE, tck = tck,
       mgp = ymgp, adj = 1, las = 1)
  axis(side = 4, at = yAt, labels = if(labelRightTicks) yTickLabels else FALSE,
       lwd = axisWidth, cex.axis = yCex, tick = TRUE, tck = tck,
       mgp = ymgp, adj = 0, las = 1)

  ## Extend the top ticks to partially enclose the plot
  if(extendTopTick){
    toptck <- ifelse(tck, 2 * tck, 0.05)
    topAt <- par()$usr[4]
    axis(side = 2, at = topAt, labels = FALSE, tick = TRUE, lwd = axisWidth, tck = toptck)
    axis(side = 4, at = topAt, labels = FALSE, tick = TRUE, lwd = axisWidth, tck = toptck)
  }

  ## labels at top of y axes
  if(nLeftTopLabels)
    for(i in 1:nLeftTopLabels)
      mtext(leftTopLabel[i], side = 3, adj = topLabelAlign,
            cex = labCex, at = xRange[1],
            line = (nLeftTopLabels - i) * labCex + 0.1)
  if(nRightTopLabels)
    for(i in 1:nRightTopLabels)
      mtext(rightTopLabel[i], side = 3, adj = 1 - topLabelAlign,
            cex = labCex, at = xRange[2],
            line = (nRightTopLabels - i) * labCex + 0.1)

  ## Inside labels at top of y axes
  offset <- innerOffset * diff(xRange)
  if(nLeftInsideLabels)
    for(i in 1:nLeftInsideLabels)
      mtext(leftInsideLabel[i], side = 3, adj = 0,
            cex = innerCex, at = xRange[1] + offset,
            line = -((1 + i) * innerCex + 0.1))
  if(nRightInsideLabels)
    for(i in 1:nRightInsideLabels)
      mtext(rightInsideLabel[i], side = 3, adj = 1,
            cex = innerCex, at = xRange[2] - offset,
            line = -((1 + i) * innerCex + 0.1))

  ## legend placement defaults
  legend.x <- 0.9 * xRange[1] + 0.1 * xRange[2]
  legend.y <- (par("usr")[4]) -
    (par("cxy")[2]) * innerCex * (2 + nLeftInsideLabels)
  
  ## labels for x axis
  if(nXAxisLabels)
    for(i in 1:nXAxisLabels)
      mtext(xAxisLabel[i], side = 1,
            line = xSpace - 1 + (i - nXAxisLabels)*labCex,
            adj = 0.5, cex = labCex)

  ## label for left axis
  if(nLeftAxisLabels)
    for(i in 1:nLeftAxisLabels)
      mtext(leftAxisLabel[i], side = 2,
            line = leftSpace - (1 + nLeftAxisLabels - i)*labCex,
            adj = 0.5, cex = labCex)

  ## label for right axis
  if(nRightAxisLabels)
    for(i in 1:nRightAxisLabels)
      mtext(rightAxisLabel[i], side = 4,
            line = rightSpace - (1 + nRightAxisLabels - i)*labCex,
            adj = 0.5, cex = labCex)

  z <- list(xy = xy, plotType = plotType, lineType = lineType,
            color = color, plotChar = plotChar, lineWidth = lineWidth,
            x = legend.x, y = legend.y, xRange = xRange, yRange = yRange,
            innerCex = innerCex, par = par())

  ## restore par() parameters except for "plt" which must keep its
  ## new value
  par(append(par.old, list(plt = par("plt"))))
  class(z) <- "scatterPlot"
  assign("latestPlot", value = z, envir = parent.frame())
  invisible(z)
}

