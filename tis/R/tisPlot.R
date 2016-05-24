tisPlot <- function(## series args
                    ...,
                    leftAxis = TRUE,
                    plotType = "l",
                    lineType = "solid",
                    lineWidth = 1.5,
                    plotChar = "*",
                    dataCex = 1,
                    color = 1,
                    midPoints = TRUE,
                    dropNA = FALSE,
                    xOffset = 0,
                    ## x axis
                    xAxisMin = NULL,
                    xAxisMax = NULL,
                    xExpandBy = 0.04,
                    xTickFreq = "Auto",
                    xTickSkip = 0,
                    xUnlabeledTickFreq = "None",
                    xUnlabeledTickSkip = 0,
                    xMinorTickFreq = "None",
                    xMinorTickSkip = 0,
                    dateFormat = "Auto",
                    xCex = 1,
                    midLabels = FALSE,
                    yearLabels = FALSE,
                    xAxisLabel = NULL,
                    xSpace = 4,
                    ## y axis
                    log = FALSE,
                    leftAxisMin = NULL,
                    leftAxisMax = NULL,
                    leftExpandBy = 0.04,
                    leftTicks = 5,
                    leftTickLocations = NULL,
                    labelLeftTicks = FALSE,
                    leftTickLabels = NULL,
                    rightAxisMin = NULL,
                    rightAxisMax = NULL,
                    rightExpandBy = 0.04,
                    rightTicks = 5,
                    rightTickLocations = NULL,
                    labelRightTicks = TRUE,
                    rightTickLabels = NULL,
                    yCex = 1,
                    extendTopTick = TRUE,
                    ## text items
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
                    innerLine = 0.5,
                    innerOffset = 0.05,
                    innerCex = 0.8,
                    foot = NULL,
                    footColor = "black",
                    footAlign = 0,
                    footCex = 0.8,
                    footSpace = -1,
                    ## layout
                    tck = 0.03,
                    axisWidth = 2,
                    start = 0,
                    end = 0,
                    boxType = "u",
                    leftMargin = -1,
                    rightMargin = -1,
                    topMargin = -1,
                    bottomMargin = -1,
                    ## shade date ranges polygons
                    nberShade = FALSE,
                    shadeDates = NULL, 
                    shadeColor = "gray",
                    shadeBorder = FALSE,
                    polyArgs = list()){
  ## Reset some graphical parameters
  parOnEntry <- par(xaxs = "i", yaxs = "i", mex = 1)
  
  if(par("new") == FALSE) plot.new()
  
  ## process (...) args into a list of named univariate series
  args <- list(...)
  serNames <- as.character(match.call(expand.dots = FALSE)$...)
  
  series <- list()
  argNames <- character(0)
  for(i in seq(along = args)){
    arg <- args[[i]]
    if(!(is.tis(arg) || is.ts(arg))){
      argName <- names(args)[i]
      msg <- paste("argument", if(argName != "") argName else serNames[i],
                   "is not a time series and will be skipped")
      warning(msg, immediate = TRUE)
      next
    }
    args.i <- as.tis(arg)
    arg.counter <- length(series)
    if(is.matrix(args.i)){
      i.cols <- dim(args.i)[2]
      for(j in 1:i.cols){
        arg.counter <- arg.counter + 1
        series[[arg.counter]] <- args.i[, j]
      }
      if(is.null(i.name <- dimnames(args.i)[[2]]))
        i.name <- paste(serNames[i], "[,", 1:i.cols, "]", sep = "")
    }
    else{
      arg.counter <- arg.counter + 1
      series[[i]] <- args.i
      i.name <- as.character(serNames[i])
    }
    
    argNames <- c(argNames, i.name)
  }
  nSeries <- length(series)
  if(nSeries == 0) stop("No valid series to plot")
  names(series) <- argNames
  
  ## Default values for some args are actually computed on the fly.  They
  ## could have been specified that way in the argument list, but that would
  ## make creating a user interface to "drive" the function more difficult.
  ## They are set here instead.
  if(missing(cex))      cex      <- par("cex")
  if(missing(lineType)) lineType <- 1:nSeries

  ## figure which series go on left and right axis
  leftAxis <- rep(leftAxis, length.out = nSeries)
  rightAxis <- !leftAxis

  ## use window if start or end is sent as argument
  if(length(start) == 1 && start == 0) start <- NULL
  if(length(end)   == 1 && end   == 0) end   <- NULL
  if((!is.null(start)) || (!is.null(end))){
    series <- lapply(series, "window", start = start, end = end, noWarn = TRUE)
  }
  ## keep starts, ends, times, and lengths
  ## for future reference
  lengths  <- sapply(series, length)
  starts   <- sapply(series, start)
  ends     <- sapply(series, end)
  tifNames <- sapply(series, tifName)
  offsets  <- rep(ifelse(midPoints, 0.5, 1), length.out = nSeries)
  xOffset  <- rep(xOffset, length.out = nSeries)
  xcts     <- vector("list", nSeries)
  for(i in 1:nSeries)
    xcts[[i]] <- POSIXct(ti(series[[i]]), offset = offsets[i]) + xOffset[i]
  
  ## Scale specific cex's by overall cex
  dataCex  <- cex * dataCex
  headCex  <- cex * headCex
  subCex   <- cex * subCex
  labCex   <- cex * labCex
  footCex  <- cex * footCex
  xCex     <- cex * xCex
  yCex     <- cex * yCex
  innerCex <- cex * innerCex

  dataCex   <- rep(dataCex,   length.out = nSeries)
  lineType  <- rep(lineType,  length.out = nSeries)
  plotType  <- rep(plotType,  length.out = nSeries)
  lineWidth <- rep(lineWidth, length.out = nSeries)
  color     <- rep(color,     length.out = nSeries)
  plotChar  <- rep(plotChar,  length.out = nSeries)
  dropNA    <- rep(dropNA,    length.out = nSeries)
  
  ## Plotting is done in the following order:
  ##  the figure margins (axes, titles, footnotes)
  ##  the lines and/or points representing the data
  ##  legends
  
  expandRange <- function(inRange, expandBy = 0.04, log = FALSE){
    ## local function
    z <- unclass(inRange)
    if(length(z) != 2) stop("length(inRange) != 2")
    if(!log){
      addFactor <- expandBy * (z[2] - z[1])
      ans <- z + c(-addFactor, addFactor)
    }
    else
      ans <- exp(Recall(log(z), expandBy, log = FALSE))
    class(ans) <- class(inRange)
    return(ans)
  }

  ## set xRange, leftRange, and rightRange
  xDataRange <- range(unlist(xcts))
  xRange <- expandRange(xDataRange, xExpandBy)
  class(xRange) <- class(xDataRange) <- c("POSIXct", "POSIXt")
  if(!is.null(xAxisMin))
    xRange[1] <- POSIXct(ti(xAxisMin, tif = "daily"), offset = 0)
  if(!is.null(xAxisMax))
    xRange[2] <- POSIXct(ti(xAxisMax, tif = "daily"), offset = 1)
  
  if(any(leftAxis)){
    leftDataRange <- range(unlist(series[leftAxis]), na.rm = TRUE)
    leftRange <- expandRange(leftDataRange, leftExpandBy, log)
    if(!is.null(leftAxisMin)) leftRange[1] <- leftAxisMin
    if(!is.null(leftAxisMax)) leftRange[2] <- leftAxisMax
  }
  if(any(rightAxis)){
    rightDataRange <- range(unlist(series[rightAxis]), na.rm = TRUE)
    rightRange <- expandRange(rightDataRange, rightExpandBy, log)
    if(!is.null(rightAxisMin)) rightRange[1] <- rightAxisMin
    if(!is.null(rightAxisMax)) rightRange[2] <- rightAxisMax
  }
  
  if(any(leftAxis) && !any(rightAxis))
    rightRange <- leftRange
  if(any(rightAxis) && !any(leftAxis))
    leftRange <- rightRange

  leftAt <- rightAt <- NULL
  ## figure out leftAt based on leftTicks and leftTickLocations
  if(!is.null(leftTickLocations))
    leftAt <- evalOrEcho(leftTickLocations)
  if(is.null(leftAt)){
    if(leftTicks > 0){
      leftAt <- signif(pretty(leftRange, n = leftTicks))
      leftAt <- leftAt[leftAt >= leftRange[1]]
    }
  }
  ## figure out rightAt based on rightTicks and rightTickLocations
  if(!is.null(rightTickLocations))
    rightAt <- evalOrEcho(rightTickLocations)
  if(is.null(rightAt)){
    if(rightTicks > 0){
      rightAt <- signif(pretty(rightRange, n = rightTicks))
      rightAt <- rightAt[rightAt >= rightRange[1]]
    }
  }
  
  if(any(leftAxis) && !any(rightAxis) && is.null(rightAt))
    rightAt <- leftAt
  if(any(rightAxis) && !any(leftAxis) && is.null(leftAt))
    leftAt <- rightAt

  if(extendTopTick){ ## extend leftRange and rightRange to include max tick locations
    if(!is.null(leftAt))   leftRange[2] <- max(c(leftRange[2], leftAt))
    if(!is.null(rightAt)) rightRange[2] <- max(c(rightRange[2], rightAt))
  }

  ## Set tick labels
  if(labelLeftTicks){
    if(is.null(leftTickLabels))
      leftTickLabels <- format(signif(leftAt, 5))
  }
  else leftTickLabels <- FALSE

  if(labelRightTicks){
    if(is.null(rightTickLabels))
      rightTickLabels <- format(signif(rightAt, 5))
  }
  else rightTickLabels <- FALSE
  
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

  nHeads             <- length(head)
  nSubs              <- length(sub)
  nFoots             <- length(foot)
  nLeftTopLabels     <- length(leftTopLabel)
  nRightTopLabels    <- length(rightTopLabel)
  nLeftInsideLabels  <- length(leftInsideLabel)
  nRightInsideLabels <- length(rightInsideLabel)
  nXAxisLabels       <- length(xAxisLabel)
  
  ylabSpace <- labCex * max(c(nLeftTopLabels, nRightTopLabels))
  headSpace <- headCex * nHeads
  subSpace <- subCex * nSubs
  topSpace <- max(c(ylabSpace, headSpace + subSpace))
  
  ## The bottom figure margin gets 4 lines of space for ticklabels,
  ## yearLabels, and/or xAxisLabel, plus whatever space is needed for
  ## footnotes.
  if(footSpace < 0) footSpace <- nFoots
  bottomSpace <- xSpace * labCex + footSpace * footCex

  csi <- par("csi")
  mex <- par("mex")
  if(is.character(leftTickLabels))
      leftTickLabelSpace <- max(strwidth(leftTickLabels, units = "inches",
                                         cex = yCex)/(csi*mex))
  else
    leftTickLabelSpace <- 0

  if(is.character(rightTickLabels))
      rightTickLabelSpace <- max(strwidth(rightTickLabels, units = "inches",
                                          cex = yCex)/(csi*mex))
  else
    rightTickLabelSpace <- 0

  mar <- par("mar")
  if( leftMargin  < 0)   leftMargin <- max(1 +  leftTickLabelSpace, mar[2])
  if(rightMargin  < 0)  rightMargin <- max(1 + rightTickLabelSpace, mar[4])
  if(topMargin    < 0)    topMargin <- max(topSpace, 1)
  if(bottomMargin < 0) bottomMargin <- bottomSpace
  
  par(mar = c(bottomMargin, leftMargin, topMargin, rightMargin))

  ## Now that margins have been set, set the plot window's user coords
  logArg <- ifelse(log, "y", "")
  plotWindow(time(xRange), leftRange, log = logArg)

  tickLocations <- function(freqString, skip = 0){
    ## local function to figure out x-axis tick locations
    ## Note: if freqString == "auto", skip is calculated by the function,
    ## ignoring the value (if any) that was passed in.
    tif <- tolower(freqString)
    if(tif == "none") return(FALSE)
    if(tif == "auto"){
      skip <- 0
      tif <- tifNames[1]
      dfreq <- tif2freq(tif)
      d <- unclass(xRange[2]) - unclass(xRange[1])
      
      if(dfreq <= 2 || d > 6*365*86400){
        locations <- POSIXct(jul(pretty(time(xRange)))-1/24)
        return(locations[between(locations, xRange[1], xRange[2])])
      }
      if(dfreq <= 4 || d > 4*365*86400){
        ## 4 to 6 year range -- tick every other quarter
        tif <- "quarterly"
        skip <- 1
      }
      else if(dfreq <= 4  || d > 365*86400) tif <-  "quarterly"
      else if(dfreq <= 12 || d > 123*86400) tif <- "monthly"
      
      if(d < 8*86400){
        if(     d <    10)   tif <- secondly(1)
        else if(d <    30)   tif <- secondly(5)
        else if(d <    60)   tif <- secondly(10)
        else if(d <   120)   tif <- secondly(15)
        else if(d <   240)   tif <- secondly(30)
        else if(d <   480)   tif <- minutely(1)
        else if(d <  1800)   tif <- minutely(5)
        else if(d <  3600)   tif <- minutely(10)
        else if(d <  7200)   tif <- minutely(15)
        else if(d < 14400)   tif <- minutely(30)
        else if(d < 21600)   tif <- hourly(1)
        else if(d < 43200)   tif <- hourly(2)
        else if(d < 86400)   tif <- hourly(4)
        else if(d < 3*86400) tif <- hourly(12)
        else                 tif <- "daily"
      }
    }
    tiRange <- c(ti(xRange[1] - 1, tif), ti(xRange[2] + 1, tif))
    if(skip > 0){
      skp1 <- skip + 1
      if(tif == "monthly"){
        tiRange[1] <- tiRange[1] + skp1 - (month(tiRange[1]) %% skp1)
        tiRange[2] <- tiRange[2] - (month(tiRange[2]) %% skp1)
      }
      if(tif == "quarterly"){
        tiRange[1] <- tiRange[1] + skp1 - (quarter(tiRange[1]) %% skp1)
        tiRange[2] <- tiRange[2] - (quarter(tiRange[2]) %% skp1)
      }
    }
    locations <- POSIXct(seq(from = tiRange[1], to = tiRange[2], by = skip + 1))
    locations[between(locations, xRange[1], xRange[2])]
    ## end of internal function
  }
  
  xTicks          <- tickLocations(xTickFreq,          xTickSkip)
  xUnlabeledTicks <- tickLocations(xUnlabeledTickFreq, xUnlabeledTickSkip)
  xMinorTicks     <- tickLocations(xMinorTickFreq,     xMinorTickSkip)

  if(length(xTicks) == 1 && as.numeric(xTicks) == 0)
    xTickLabels <- FALSE
  else {
    if(tolower(dateFormat) == "auto"){
      d <- unclass(xRange[2]) - unclass(xRange[1])
      if(     d <      1.1 * 60)           dateFormat <- "%S"
      else if(d <      1.1 * 60 * 60)      dateFormat <- "%M:%S"
      else if(d <      1.1 * 60 * 60 * 24) dateFormat <- "%H:%M"
      else if(d <      2   * 60 * 60 * 24) dateFormat <- "%a %H:%M"
      else if(d <      7   * 60 * 60 * 24) dateFormat <- "%a"
      else if(d <     50   * 60 * 60 * 24) dateFormat <- "%b %d"
      else if(d < 365 *1.1 * 60 * 60 * 24) dateFormat <- "%b"
      else                                 dateFormat <- "%Y"
    }
    xTickLabels <- format(jul(xTicks), format = dateFormat)
  }

  ## Place titles and footnotes
  ##   head
  if(nHeads)
    for(i in 1:nHeads)
      mtext(head[i], side = 3, adj = headAlign, cex = headCex,
            line = headSpace + subSpace - i * headCex, font = par("font.main"))
  ##   sub
  if(nSubs)
    for(i in 1:nSubs)
      mtext(sub[i], side = 3, adj = headAlign, cex = subCex,
            line = subSpace - i * subCex, font = par("font.sub"))
  ##   foot
  if(nFoots > 0){
    footColor <- rep(footColor, length.out = nFoots)
    for(i in 1:nFoots)
      mtext(foot[i], side = 1, adj = footAlign, cex = footCex,
            line = xSpace * labCex - 1 + i * footCex, col = footColor[i])
  }
  ## shaded date ranges
  if(nberShade && is.null(shadeDates)) shadeDates <- nberDates()
  if(!is.null(shadeDates)){
    ymdShadeArgs <- list(ymds = ymd(shadeDates),
                         col = shadeColor,
                         border = shadeBorder,
                         xrange = time(xDataRange))
    if(length(polyArgs) > 0)
      ymdShadeArgs <- updateList(ymdShadeArgs, polyArgs)
    do.call("ymdShade", ymdShadeArgs)
  }

  ## draw the X axis
  box(lwd = axisWidth, bty = boxType)

  ## mgp parameter is supposed to be the margin line in mex units for
  ## c(axis title, axis labels, axisline).  We don't use the first, and the
  ## last is always 0, but we need to compute a reasonable number for the
  ## second.
  xmgp <- par("mgp")
  xmgp[2] <- max(0.6*xCex - 0.5, 0)
  
  if(!(length(xUnlabeledTicks) == 1 && as.logical(xUnlabeledTicks) == FALSE))
    axis(side = 1, at = time(xUnlabeledTicks), labels = FALSE, tick = TRUE,
         lwd = axisWidth, tck = tck, mgp = xmgp)
  
  if(!(length(xMinorTicks) == 1 && as.logical(xMinorTicks) == FALSE))
    axis(side = 1, at = time(xMinorTicks), labels = FALSE, tick = TRUE,
         lwd = 1, tck = tck*2/3, mgp = xmgp)
  
  if(!(length(xTicks) == 1 && as.logical(xTicks) == FALSE)){
    if(midLabels){
      ## do the ticks, then the labels
      axis(side = 1, at = time(xTicks), labels = FALSE, tick = TRUE,
           lwd = axisWidth, tck = tck, mgp = xmgp)
      
      if(length(xTicks) > 1){
        ucxt <- unclass(xTicks)
        labat <- structure((ucxt[-1] + ucxt[-length(xTicks)])/2,
                           class = c("POSIXt", "POSIXct"))
        axis(side = 1, at = time(labat), labels = xTickLabels[-1],
             tick = FALSE, mgp = xmgp, cex.axis = xCex)
      }
      else {
        labat <- xTicks
        axis(side = 1, at = time(labat), labels = xTickLabels,
             tick = FALSE, mgp = xmgp, cex.axis = xCex)
      }
    }
    else {
      if(midPoints[1]){
        axis(side = 1, at = time(xTicks), labels = FALSE, tick = TRUE, lwd = axisWidth,
             tck = tck, mgp = xmgp)
        mtext(text = xTickLabels, side = 1, line = xmgp[2], at = time(xTicks),
              adj = 0.9, cex = xCex)
      }
      else 
        axis(side = 1, at = time(xTicks), labels = xTickLabels, tick = TRUE,
             lwd = axisWidth, tck = tck, mgp = xmgp, cex.axis = xCex)      
    }
  }
 ## draw yearLabels on line below xTickLabels
  if(yearLabels){
    yt <- sort(c(xDataRange, tickLocations("annual", skip = 0)))
    nyt <- length(yt)
    yearAt <- yt[-nyt] + unclass(difftime(yt[-1], yt[-nyt], units = "secs")/2)
    yearLabels <- format(yearAt, "%Y")
    if(length(xTickLabels) == 1 && !is.character(xTickLabels) && !as.logical(xTickLabels))
      yearLabelLine <- 0
    else
      yearLabelLine <- labCex
    axis(side = 1, at = time(yearAt), labels = yearLabels, tick = FALSE,
         line = yearLabelLine, mgp = xmgp, cex.axis = xCex)
  }

  ## drop NA observations from the plotting lists
  px <- lapply(xcts, time)
  py <- series
  if(any(dropNA)){
    for(i in(1:length(px))){
      if(dropNA[[i]]){
        naSpots <- is.na(py[[i]])
        px[[i]] <- px[[i]][!naSpots]
        py[[i]] <- py[[i]][!naSpots]
      }
    }
  }
  
  ## plot the leftscale timeseries and axis
  if(any(leftAxis))
    for(i in(1:nSeries)[leftAxis])
      lines(px[[i]], py[[i]], 
            type = plotType[i], lty = lineType[i], lwd = lineWidth[i],
            col = color[i], pch = plotChar[i], cex = dataCex[i])

  if(!(length(leftAt) == 1 && as.logical(leftAt) == F))
    axis(side = 2, at = leftAt, labels = leftTickLabels, lwd = axisWidth,
         tick = TRUE, tck = tck, mgp = c(3, 0.5, 0), adj = 1, las = 1, cex.axis = yCex)
  ## plot right scale timeseries and axis
  if(any(rightAxis)){
    plotWindow(time(xRange), rightRange, log = logArg)
    for(i in(1:nSeries)[rightAxis])
      lines(px[[i]], py[[i]], 
            type = plotType[i], lty = lineType[i], lwd = lineWidth[i],
            col = color[i], pch = plotChar[i], cex = dataCex[i])
  }
  if(!(length(rightAt) == 1 && as.logical(rightAt) == F))
    axis(side = 4, at = rightAt, labels = rightTickLabels, lwd = axisWidth,
         tick = TRUE, tck = tck,  mgp = c(3, 0.5, 0), adj = 0, las = 1, cex.axis = yCex)
  ## Extend the top ticks to partially enclose the plot
  if(extendTopTick){
    toptck <- ifelse(tck, 2 * tck, 0.05)
    topAt <- par()$usr[4]
    if(par("ylog")) topAt <- 10^topAt
    axis(side = 2, at = topAt, labels = FALSE, tick = TRUE, lwd = axisWidth, tck = toptck)
    axis(side = 4, at = topAt, labels = FALSE, tick = TRUE, lwd = axisWidth, tck = toptck)
  }
  ## labels at top of y axes
  if(nLeftTopLabels)
    for(i in 1:nLeftTopLabels)
      mtext(leftTopLabel[i], side = 3, adj = topLabelAlign,
            cex = labCex, at = time(xRange[1]),
            line = (nLeftTopLabels - i) * labCex + 0.1)
  if(nRightTopLabels)
    for(i in 1:nRightTopLabels)
      mtext(rightTopLabel[i], side = 3, adj = 1 - topLabelAlign,
            cex = labCex, at = time(xRange[2]),
            line = (nRightTopLabels - i) * labCex + 0.1)
  ## Inside labels at top of y axes
  offset <- innerOffset * diff(time(xRange))
  if(nLeftInsideLabels)
    for(i in 1:nLeftInsideLabels)
      mtext(leftInsideLabel[i], side = 3, adj = 0,
            cex = innerCex, at = time(xRange[1]) + offset,
            line = -((innerLine + i) * innerCex + 0.1))
  if(nRightInsideLabels)
    for(i in 1:nRightInsideLabels)
      mtext(rightInsideLabel[i], side = 3, adj = 1,
            cex = innerCex, at = time(xRange[2]) - offset,
            line = -((innerLine + i) * innerCex + 0.1))
  ## yLegendOffset is in relative ([0..1]) units
  yLegendOffset <- (par("cin")/par("pin"))[2] * innerCex * (2 + nLeftInsideLabels)

  ## label for x axis
  if(nXAxisLabels)
    for(i in 1:nXAxisLabels)
      mtext(xAxisLabel[i], side = 1,
            line = xSpace - 1 + (i - nXAxisLabels)*labCex,
            adj = 0.5, cex = labCex)
  
  z <- list(series = series, dateFormat = dateFormat,
            plotType = plotType, lineType = lineType, dataCex = dataCex,
            color = color, plotChar = plotChar, lineWidth = lineWidth,
            yLegendOffset = yLegendOffset, cex = cex,
            xRange = time(xRange), leftRange = leftRange, rightRange = rightRange,
            midPoints = midPoints, par = par())

  ## restore par() parameters except for "plt" which must keep its 
  ## new value
  par(append(parOnEntry, list(plt = par("plt"))))
  class(z) <- "tisPlot"
  assign("latestPlot", value = z, envir = parent.frame())
  invisible(z)
}

tisLegend <- function(..., xrel = 0.1, yrel = 0.1, xjust = 0, yjust = 1,
                      boxType = "n", ncol = 1, cex = 1){
  ## assumes the existence of latestPlot, a tisPlot.
  if(exists("latestPlot", envir = parent.frame()))
    sp <- get("latestPlot", envir = parent.frame())
  else {
    if(exists("latestPlot", envir = globalenv()))
      sp <- get("latestPlot")
    else
      stop("No latestPlot to label")
  }
  
  legendArgs <- list(legend = names(sp$series),
                     col = sp$color,
                     lty = sp$lineType,
                     lwd = sp$lineWidth,
                     pch = sp$plotChar,
                     bty = boxType,
                     cex = cex*sp$cex,
                     ncol = ncol,
                     xjust = xjust,
                     yjust = yjust)
  ## handle the other arguments, if there are any
  arglist <- list(...)
  args.n <- length(arglist)
  if(args.n){
    ## if user called " tisLegend(locator(), ...) ", arglist will
    ## have some list components, which have to be extracted
    list.args <- which(sapply(arglist, is.list))
    if(length(list.args) > 0)
      arglist <- c(arglist[[list.args]], arglist[-list.args])
    argnames <- names(arglist)
    ## Update legendArgs with the components of arglist, replacing where
    ## appropriate
    legendArgs[argnames] <- arglist
  }
  usr <- par("usr")
  if(is.character(arglist$x))
    legendArgs$x <- arglist$x
  else {
    if(is.null(arglist$x)){
      legendArgs$x <- usr[2]*xrel + usr[1]*(1 - xrel)
      if(par("xlog") && !is.character(legendArgs$x)) legendArgs$x <- 10^legendArgs$x
    }
    else legendArgs$x <- arglist$x
    
    if(is.null(arglist$y)){
      if(missing(yrel)) yrel <- sp$yLegendOffset
      legendArgs$y <- usr[4]*(1 - yrel) + usr[3]*yrel
      if(par("ylog")) legendArgs$y <- 10^legendArgs$y
    }
    else legendArgs$y <- arglist$y
  }
  ## Make certain args the same length as legend
  n <- length(legendArgs$legend)
  legendArgs$col <- rep(legendArgs$col, length.out = n)
  legendArgs$lty <- rep(legendArgs$lty, length.out = n)
  legendArgs$lwd <- rep(legendArgs$lwd, length.out = n)
  legendArgs$pch <- rep(legendArgs$pch, length.out = n)
  type <- rep(sp$plotType, length.out = n)
  hasPoints <- (type == "p" | type == "b" | type == "o")
  hasLines  <- (type == "l" | type == "b" | type == "c" |
                type == "o" | type == "h")
  if(all(! hasPoints ))
    legendArgs$pch <- NULL
  else
    legendArgs$pch[!hasPoints] <- ""

  if(is.numeric(legendArgs$lty))
    legendArgs$lty[!hasLines] <- 0
  else
    legendArgs$lty[!hasLines] <- "blank"
  if(all(! hasLines)) legendArgs$lwd <- NULL

  ## color the legend text
  if(is.null(legendArgs$text.col))
    legendArgs$text.col <- legendArgs$col
  legendArgs$text.col <- rep(legendArgs$text.col, length.out = n)
  
  do.call("legend", legendArgs)
  class(legendArgs) <- "tisLegend"
  invisible(legendArgs)
}
