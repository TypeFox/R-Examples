## $Id: utilities.R 6053 2015-02-22 20:23:45Z bpikouni $ 

## Functions & objects internal for the cg Library
## and definitely not intended for user-level calls

#####################

stndErr <- function(x) {
  ##
  ## PURPOSE: Calculate standard error measure from a group of data.
  ## NOTE: We use of "n" in the function name
  ## instead of "stderr" because SV4 "stderr"
  ## is a system base function that is important in connections.
  x <- x[!is.na(x)]
  sd(x) / sqrt(length(x))
}

#####################

geoMean <- function(x, na.rm=TRUE) {
  ##
  ## PURPOSE: Calculate the geometric mean from a group of data
  if(na.rm) {
    x <- x[!is.na(x)]
  }
  if(any(x <= 0)) {
    stop(cgMessage("There appear to be some non-positive",
                   "values in the data. All values",
                   "must be greater than zero."))
  }
  exp(mean(log(x)))
}

#####################

pctToRatio <- function(x) { (x/100) + 1 }

#####################

makeZeroScore <- function(x) {
  respvec <- unique(stripmiss(x[x > 0]))
  return(exp(spline(respvec, log(respvec), method="natural",
                    xmin=0, n=length(respvec) + 1)$y[1]))
}

#####################

unwind <- function(x) {
  as.vector(unlist(x))
}

unwrap <- unwind

#####################

graphStampCG <- function(grid=TRUE) {
  ##
  ## PURPOSE: ID Stamp on graphs
  ##
  validBoolean(grid)
  ## Reads off package object "cgVersion"
  dateitems <- unlist(strsplit(date(), split=" "))
  ## Remove extraneous string if present due to single digit days
  dateitems <- dateitems[dateitems!=""]

  ## sequence: daynumber, month, year, time, weekday
  newdate <- paste(dateitems[3]," ",dateitems[2]," ",
                   dateitems[5],"  \n",dateitems[4]," local  \n",sep="")
  thetext <- paste("cg",
            " v", cgVersion, "  \n",
            newdate,"", sep="")
  if(grid) {
    grid::grid.text(thetext, x = 0.99, y = 0.01, 
                    default.units = "npc", 
                    just = c("right", "bottom"),
                    gp=gpar(cex=0.6))
  }
  else { ## not grid
    mtext(text=thetext,
          side = 1,
          line = -1, outer = TRUE,
          adj = 1, cex = 0.6)
  }
  invisible()
}

#####################

setupAxisTicks <- function(x, ratio=FALSE, percent=FALSE,
                           difference=FALSE, axis="y",
                           logscale=TRUE, digits=NULL, remticks=TRUE,
                           offset=NULL, grid=FALSE, xcex=1, ycex=1) {
  ##
  ## PURPOSE:  Try to setup good ticks for log spacing
  ## or original scale spacing: 
  ## Ensure no overlap, and
  ## Ensure space for later addition of min and max
  ## NOTE: x is assumed to span the whole axis range of data
  ##
  if(is.null(digits)) digits <- 4
  
  x <- x[!is.na(x)]

  if (!is.null(offset)) {
    if(logscale==FALSE) {
      stop(cgMessage("The offset argument must be set to NULL when",
                     "logscale=FALSE, as they are not compatible."))
    }
    x <- (x - offset)
  }
  
  if(any(x < 0)) {
    ## generating ticks on the magnitude of values
    x.pos <- abs(x)
  }

  else {
    x.pos <- x
  }
  
  log10endpt <- log10(x.pos) 

  if(ratio==TRUE) {
    ## to impose that zero change is always represented
    log10endpt <- c(0, log10endpt)
  }

  digitplacement <- floor(max(log10endpt))
  log10ticks <- pretty(log10endpt)
  tickmarks <- round(10^log10ticks, -digitplacement)

  ## need to filter extreme values for case where logscaling
  ## essentially has no effect, i.e. the range is clearly small.
  tickmarks <- tickmarks[tickmarks >= round(min(x), digits) &
                         tickmarks <= round(max(x), digits)]

  ## Also ensure that powers of 10 from 1 thru range are plotted
  powersof10 <- 10^(0:digitplacement)

  ## Additional adjustments for Ratio/Percent Change
  if(ratio==TRUE) {
    powersof10 <- (powersof10/10) + 1  
  }
  
  powersof10 <- powersof10[powersof10 >= round(min(x), digits) &
                           powersof10 <= round(max(x), digits)]

  ## and that some interim values that are familiar are there
  if(ratio==TRUE) {
    interimticks <- as.vector(outer(c(-9.9,-9,-7.5,-5,-2.5,-1,1,2.5,5),
                                    10^(0:digitplacement)))
    ## Additional adjustments for percent change to ratio
    interimticks <- (interimticks/10) + 1
  }
  else {
    interimticks <- as.vector(outer(c(1,2,5),10^(0:digitplacement)))    
  }
  
  ## use only those that fall within range
  interimticks <- interimticks[interimticks >= (round(min(x), digits) - 1e-8) &
                               interimticks <= (round(max(x), digits) + 1e-8) ]
  
  ## If Ratio, we also need to include the
  ## "no-difference" reference which is 1 on the ratio scale
  ## Also need to make sure extremes of axes are included
  if(ratio==TRUE) {
    ## Cover the two extremes of the axes via percent scaling because we can
    ## work with integers
    ends <- c(floor(100*(min(x.pos) - 1)),
              round(100*(min(x.pos) - 1),0),
              round(100*(max(x.pos) - 1),0),
              ceiling(100*(max(x.pos) - 1)))
    
    ## if 0 at lower end replace with 1
    ends[ends==0] <- 1

    ## transform to ratio
    ends <- (ends/100) + 1
    
    xticks <- c(1, tickmarks, powersof10, interimticks, ends,
                pretty(x.pos))
  }
  else {
    ## Cover the two extremes of the axes
    ends <- round(range(x), digits)

    xticks <- c(tickmarks, powersof10, interimticks, ends,
                pretty(x))
    ## pretty(x.pos))
  }

  ## Remove 0 if it has crept in there if we plan use for log-spacing
  if(logscale==TRUE) {
    xticks <- xticks[xticks > 0]
  }

  ## Annotation in case the original data had negative values
  ## Except do not reflect the min and max.
  if(any(x < 0)) {
    xticks <- c(-xticks, 0, xticks)
    xticks <- xticks[!(xticks==round(-min(x), digits) | xticks==round(-max(x),
                                                          digits))]
  }

  ## If ratio, then try to ensure symmetry of tickmarks in the range of 0.1 to
  ## 1.9
  if(ratio==TRUE) {
    upperbound <- max(min(round(x.pos,1), 1.9), 1.1)
    lowerbound <- min(min(round(x.pos,1)), 0.9)
    if(upperbound < 1.8) upperbound <- upperbound + 0.2
    if(lowerbound > 0.2) lowerbound <- lowerbound - 0.2
    xticks <- c(seq(lowerbound, upperbound, by=0.1), xticks)
  }

  ## If difference, then try to ensure symmetry of tickmarks
  if(difference==TRUE) {
    xticks <- c(0, -xticks, xticks)
  }

  ## If ratio, then try to ensure symmetry of tickmarks for smaller ranges
  if(ratio==TRUE) {
    xticks <- c(xticks[xticks < max(xticks)] - 1, xticks)
    xticks <- xticks[xticks > 0]
  }

  ## If ratio, round off values > 2 to integers
  if(ratio==TRUE) {
    xticks[xticks > 2] <- round(xticks[xticks > 2], 0)
  }

  ## Exclude anything outside the min or max of original data, Note that 1e-8 is
  ## added as padding to insure that min and max are included in the comparison
  ## but that nothing else of consequence will be.
  xticks <- xticks[xticks >= (round(min(x), digits) - 1e-8) &
                   xticks <= (round(max(x), digits) + 1e-8)]
  
  ## Take out duplicates and ensure sorted order
  xticks <- sort(unique(xticks))
  ## Empty case handling
  if(length(xticks)==0) {
    xticks <- 1
    names(xticks) <- ""
    return(xticks)
  }

  ## Format tick labels
  ## Try to remove extraneous decimal places for larger numbers
  if(ratio==FALSE && axis=="y") {
    names(xticks) <- format(chopZeroes(round(xticks, digits)),
                            small.interval=0, justify="right")
  }
  else if(ratio==FALSE && axis=="x") {
    names(xticks) <- format(chopZeroes(round(xticks, digits)),
                            small.interval=0, justify="right")
  }
  else if(percent==TRUE && axis=="x") {
    names(xticks) <- fmtRatioToPercent(xticks)
  }
  else if(percent==TRUE && axis=="y") {
    names(xticks) <- fmtRatioToPercent(xticks)
  }
  else if(difference==TRUE && axis=="x") {
    names(xticks) <- fmtDifference(xticks)
  }
  else if(difference==TRUE && axis=="y") {
    names(xticks) <- format(chopZeroes(round(xticks, digits)), justify="right")
  }
  else if(ratio==TRUE) {
    names(xticks) <- fmtRatio(xticks)
  }

  names(xticks) <- trimWhiteSpace(names(xticks))
  ## Remove duplicate *labels* (which may have come from rounding digits)
  xticks <- xticks[which(!duplicated(names(xticks)))]
  if(!is.null(offset)) { xticks <- (xticks + offset) }

  if(remticks==TRUE) {
    ## Remove overlapping ticks
    xticks <- rmTicks(xticks, axis=axis, logscale=logscale,
                      difference=difference, percent=percent,
                      ratio=ratio,
                      digits=digits, grid=grid, xcex=xcex, ycex=ycex)
  }
  return(xticks)
}

#####################

setupLog10AxisTicks <- function(endpt, offset) {
  ##
  ## PURPOSE: Create log110 tickmarks intended for
  ## axis = 4 and that can accomodate log(x + addconstant)
  ## transformed data
  ##  
  prettyticks.log10endpt <- if(is.null(offset)) {
    pretty(log10(endpt))
  }
  else {
    pretty(log10(endpt - offset))
  }

  ticks.log10endpt <- if(is.null(offset)){
    prettyticks.log10endpt
  }
  else {
    log10(10^prettyticks.log10endpt + offset)
  }

  names(ticks.log10endpt) <- prettyticks.log10endpt
  return(ticks.log10endpt)
}


#####################

tryAgain <- paste("Please check your argument inputs or/and data, correct as needed, and",
                  "try again. Typing in traceback() next may also help.")

#####################

seeHelpFile <- function(x) {
  return(paste("The help file for", x, "may provide details."))
}


#####################

paragraphWrap <- function (x, width = 80) {
  ##
  ## PURPOSE: Format a string for use in stop() or warning().
  ## 
  ## NOTE: Adapted from R function base:::strwrap
  ##
  indent <- 0
  exdent <- 0
  prefix <- ""
  indentString <- paste(rep(" ", indent), collapse = "")
  exdentString <- paste(rep(" ", exdent), collapse = "")
  y <- list()

  z <- lapply(strsplit(x, "\n[        \n]*\n"), strsplit, "[  \n]")
  
  for (i in seq(along = z)) {
    yi <- character(0)
    for (j in seq(along = z[[i]])) {
      words <- z[[i]][[j]]
      nc <- nchar(words)
      if (any(nc == 0)) {
        zLenInd <- which(nc == 0)
        zLenInd <- zLenInd[!(zLenInd %in% (grep("\\.$", words) + 1))]
        if (length(zLenInd) > 0) {
          words <- words[-zLenInd]
          nc <- nc[-zLenInd]
        }
      }
      if (length(words) == 0) {
        yi <- c(yi, "", prefix)
        next
      }
      currentIndex <- 0
      lowerBlockIndex <- 1
      upperBlockIndex <- integer(0)
      lens <- cumsum(nc + 1)
      first <- TRUE
      maxLength <- width - nchar(prefix) - indent
      while (length(lens) > 0) {
        k <- max(sum(lens < maxLength), 1)
        if (first) {
          first <- FALSE
          maxLength <- maxLength + indent - exdent
        }
        currentIndex <- currentIndex + k
        if (nc[currentIndex] == 0) 
          upperBlockIndex <- c(upperBlockIndex, currentIndex - 
                               1)
        else upperBlockIndex <- c(upperBlockIndex, currentIndex)
        if (length(lens) > k) {
          if (nc[currentIndex + 1] == 0) {
            currentIndex <- currentIndex + 1
            k <- k + 1
          }
          lowerBlockIndex <- c(lowerBlockIndex, currentIndex + 
                               1)
        }
        if (length(lens) > k) 
          lens <- lens[-(1:k)] - lens[k]
        else lens <- NULL
      }
      nBlocks <- length(upperBlockIndex)
      s <- paste(prefix, c(indentString, rep(exdentString, 
                                             nBlocks - 1)), sep = "")
      for (k in (1:nBlocks))
        s[k] <- paste(s[k],
                      paste(words[lowerBlockIndex[k]:upperBlockIndex[k]], 
                            collapse = " "), sep = "")
      yi <- c(yi, s, prefix)
    }
    y <- c(y, list(yi[-length(yi)]))
  }
  ## We always want it " simplify'ed " (see base::strwrap)
  y <- unlist(y)

  ## Next to lastly we add the functionality of putting it into a string with
  ## linebreaks for its intended use in stop() or warning()
  y <- paste(y, collapse="\n")

  ## and add a blank space after periods followed by an upper case characters
  y <- gsub('\\.([A-Z])', '\\.\\ \\1', y, perl=TRUE)  
  
  ## Add a newline so the first line is flush left when shown
  ## and another at the end.
  return(paste("\n", y, "\n",sep=""))
}

#####################

cgMessage <- function(..., warning=FALSE) {
  paragraphWrap(paste(..., ifelse(warning==FALSE, tryAgain, "")))
}

#####################

factorInSeq <- function(x, labels) {
  ##
  ## PURPOSE: Laziness motivated wrapper for creating
  ## a factor with the order of its levels determined
  ## by the data sequence
  ##
  inseqlevels <- unique(x)
  if(missing(labels)) {
    labels <- as.character(inseqlevels)
  }
  factor(x, levels=unique(x), labels=labels)
}

#####################

setupGrpNameTicks <- function(x, grplocation=1:length(x),
                              cexinit=1, cexthreshold=0.80, grid=FALSE) {
  ##
  ## PURPOSE:  Handle group name labels so there is no overlap
  ##           Either shrink font, or rotate 45 degrees (or both)
  ##
  ## NOTE: For traditional graphics,
  ## Assume this is to be used with a graph that has just been
  ## plotted, since we need critical values from the device, which
  ## are given in inches.
  ##
  ## x assumed to be a character vector already in preferred order
  ##
  if(!grid) { ## Traditional Graphics
    axislength <- par("pin")[1]
    axismin <- par("usr")[1] 
    axisrange <- par("usr")[2] - par("usr")[1]

    ## Since widths will vary
    chrspace <-  par("cin")[1]
    chrspaces <- chrspace*nchar(x)
  }
  else { ## grid, work with npc units
    usrx <- convertX(unit(0:1, "npc"), "native", TRUE)
    usry <- convertY(unit(0:1, "npc"), "native", TRUE)
    ## grid lattice natural coordinates for panel
    axislength <- 1 

    chrspaces <- convertX(stringWidth(x), "npc", TRUE)
    chrspace <- max(convertX(stringWidth(letters), "npc", TRUE))
    axismin <- usrx[1]
    axisrange <- usrx[2] - usrx[1]

  }
  
  ## Need to map response scale to device inches scale
  n <- length(x)
  p <- (axislength/axisrange) * (grplocation - axismin)
  
  ## Construct comparisons to determine overlap
  ppart <- p[-1] - p[-n]
  chrpart <- (chrspaces[-n] + chrspaces[-1])/2
  
  comparisons <- ppart - chrpart
  
  ## Initializations
  thecex <- cexinit
  thesrt <- 0
  theadj <- 0.5
  
  overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
  while(any(overlapchecks==FALSE) && thecex >= cexthreshold) {
    thecex <- thecex - 0.05
    chrspace <- chrspace*thecex
    overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
  }

  if(thecex < cexthreshold) {
    thesrt <- 45
    theadj <- 1
    ## Now need to check overlap in the rotated scale, which
    ## can be done assuming 90 degree rotation, since
    ## character height now serves as the magic overlap msmt.
    chrspace <-
      if(!grid) {
        par("cin")[2]
      }
      else {
        max(convertX(stringHeight(letters), "inches", TRUE))
      }
    ## we only need to worry about 1 character since
    ## they labels are virtually perpendicular to the axis   
    chrspaces <- rep(chrspace, length(x))
    chrpart <- (chrspaces[-n] + chrspaces[-1])/2
    comparisons <- ppart - chrpart
    
    thecex <- 1
    overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
    while(any(overlapchecks==FALSE) && thecex >= cexthreshold) {
      thecex <- thecex - 0.05
      chrspace <- chrspace*thecex
      overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
    }
    
  }
  
  if(!grid) {
    ## If the labels are to be rotated, then we need to make sure they are
    ## shrunk within the margin available.
    ## Only worry about this for traditional graphics
    if(thesrt!=0) {
      themaxlength <- max(strwidth(x,
                          units="inches"))/sqrt(2)
      thecex <- min(thecex, (par("mai")[1]/themaxlength)*thecex)
    }
  }

  return(list(cex=thecex, srt=thesrt, adj=theadj))
  
}

xTicksCex <- function(x, thenames=names(x),
                      cexinit=1, cexthreshold=0.70,
                      grid=FALSE) {
  ##
  ##
  thenames <- names(x)
  if(!grid) { ## Traditional Graphics
    axislength <- par("pin")[1]
    axismin <- par("usr")[1] 
    axisrange <- par("usr")[2] - par("usr")[1]

    ## Since widths will vary
    chrspace <-  par("cin")[1]
    chrspaces <- chrspace*nchar(thenames)
  }
  else { ## grid lattice natural coordinates for panel
    usrx <- convertX(unit(0:1, "npc"), "native", TRUE)
    axislength <- 1
    axismin <- usrx[1]
    axisrange <- usrx[2] - usrx[1]

    chrspaces <- convertX(stringWidth(thenames), "npc", TRUE)
    chrspace <- max(convertX(stringWidth(letters), "npc", TRUE))
  }
  
  ## Need to map response scale to device inches scale for traditional,
  ## npc for grid
  n <- length(x)
  p <- (axislength/axisrange) * (x - axismin)
  
  ## Construct comparisons to determine overlap
  ppart <- p[-1] - p[-n]
  chrpart <- (chrspaces[-n] + chrspaces[-1])/2
  
  comparisons <- ppart - chrpart
  
  ## Initializations
  thecex <- cexinit
  
  overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
  while(any(overlapchecks==FALSE) && thecex > cexthreshold) {
    thecex <- thecex - 0.05
    chrspace <- chrspace*thecex
    overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
  }
  
  return(thecex)
}

yTicksCex <- function(y,
                      cexinit=1, cexthreshold=0.70,
                      grid=FALSE) {
  ##
  ##
  thenames <- names(y)
  if(!grid) { ## Traditional Graphics
    axislength <- par("pin")[2]
    axismin <- par("usr")[2] 
    axisrange <- par("usr")[4] - par("usr")[3]

    chrspace <-  par("cin")[2]
    chrspaces <- chrspace*nchar(thenames)
  }
  else {
    usry <- convertY(unit(0:1, "npc"), "native", TRUE)
    ## grid lattice natural coordinates for panel
    axislength <- 1
    axismin <- usry[1]
    axisrange <- usry[2] - usry[1]      

    chrspaces <- convertY(stringHeight(thenames), "npc", TRUE)
    chrspace <- max(convertY(stringHeight(letters), "npc", TRUE))
  }
  
  ## Need to map response scale to device inches scale
  n <- length(y)
  p <- (axislength/axisrange) * (y - axismin)
  
  ## Construct comparisons to determine overlap
  ppart <- p[-1] - p[-n]
  chrpart <- (chrspaces[-n] + chrspaces[-1])/2
  
  comparisons <- ppart - chrpart
  
  ## Initializations
  thecex <- cexinit
  
  overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
  while(any(overlapchecks==FALSE) && thecex > cexthreshold) {
    thecex <- thecex - 0.05
    chrspace <- chrspace*thecex
    overlapchecks <- c(TRUE, comparisons > chrspace/1.2)
  }
  
  return(thecex)
}

#####################


rmTicks <- function(x, axis="y", logscale=TRUE, ratio=FALSE, percent=FALSE,
                    difference=FALSE, digits=NULL, grid=FALSE, xcex=1, ycex=1) {
  ##
  ## PURPOSE:  Try to prevent overlapping tickmarks
  ##
  ## NOTE: Assume this is to be used with a graph that has just been
  ## plotted, since we need critical values from the device, which
  ## are given in inches.
  ##
  x.orig <- x
  x <- sort(x)

  if(logscale) {
    x <- log10(x)
  }
  
  if(!grid) { ## Traditional Graphics
    if(axis=="x") {
      chrspace <- par("cin")[1]
      axislength <- par("pin")[1]
      axismin <- par("usr")[1] 
      axisrange <- par("usr")[2] - par("usr")[1]
      
      chrspaces <- chrspace*nchar(names(x))
      thecex <- xcex
    }
    
    else { ## axis="y"
      chrspace <- par("cin")[2]
      axislength <- par("pin")[2]
      axismin <- par("usr")[3] 
      axisrange <- par("usr")[4] - par("usr")[3]
      
      ## in the y-axis case we only need to worry about 1 character since
      ## the y labels are perpendicular to the axis
      chrspaces <- rep(chrspace, length(x))
      thecex <- ycex
    }
  }
  else if(grid) {
    usrx <- convertX(unit(0:1, "npc"), "native", TRUE)
    usry <- convertY(unit(0:1, "npc"), "native", TRUE)
    ## grid lattice natural coordinates for panel
    axislength <- 1 
    if(axis=="x") {
      chrspaces <- convertX(stringWidth(names(x)), "npc", TRUE)
      axismin <- usrx[1]
      axisrange <- usrx[2] - usrx[1]
      chrspace <- convertX(stringWidth("0"), "npc", TRUE) 
      thecex <- xcex
    }
    
    else { ## axis="y"
      chrspaces <- convertY(stringHeight(names(x)), "npc", TRUE)
      axismin <- usry[1]
      axisrange <- usry[2] - usry[1]
      chrspace <- convertY(stringHeight("0"), "npc", TRUE) 
      thecex <- ycex
    }
  }

  ## Need to map response scale to device inches scale for traditional,
  ## npc for grid
  p <- (axislength/axisrange) * (x - axismin)
  n <- length(p)

  ## Stay in the p-space for conducting comparisons
  ## TRUE means to keep (no overlap)
  ## Keep min and max, and order statistics 2 and n - 1
  if(n > 4) {
    overlapchecks <- rep(TRUE, n)
    lastOKtickindx <- 2
    for(i in 3:(n-2)) {
      if( ( 1.0*(p[i] - p[lastOKtickindx]) -
           (thecex*(chrspaces[i] + chrspaces[lastOKtickindx])/2) ) > 0) {
        overlapchecks[i] <- TRUE
        lastOKtickindx <- i
      }
      else {
        overlapchecks[i] <- FALSE
      }
    }

    ## Make sure the order n - 1 statistic is included
    if(overlapchecks[n-1]==FALSE) {
      overlapchecks[n-2] <- FALSE
      overlapchecks[n-1] <- TRUE
    }

    ## Purposely trim the minimum and maximum
    overlapchecks[c(1, n)] <- FALSE
    
    keep <- p[overlapchecks]

    ## Go back to the x space
    out <- x[is.element(p, keep)]
  }
  else out <- x

  ## Ensure that zero/one no-difference tickmark is always included
  if(difference==TRUE || percent==TRUE)  {
    if(!any(out==0)) {
      out <- c(0, out)

      if(ratio==TRUE && percent==FALSE) {
        names(out)[1] <- "1"
      }
      else names(out)[1] <- "0"

      out <- sort(out)
    }
  }

  if(logscale) {
    out <- 10^(out)
  }

  ## If no entries
  if(all(is.na(out))) {
    out <- mean(range(x.orig))
    names(out) <- " "
  }

  return(out)

}



#####################

minmaxTicks <- function(x, theaxis="y", logscale=TRUE,
                        ratio=FALSE, percent=FALSE, difference=FALSE,
                        digits=NULL, offset=NULL, zeroscore=NULL,
                        maxonly=FALSE, grid=FALSE, cex=0.7,...) {
  ##
  ## PURPOSE:  Add min and max to axis tick labels
  ##
  if(!is.null(zeroscore) && !is.null(offset)) {
    stop(cgMessage("The zeroscore and offset arguments",
                   "cannot be specified at the same time."))
  }

  axisrange <- range(x)
  if(is.numeric(offset)) { x <- x - offset }
  labelrange <- range(x) 
  if(is.null(digits)) digits <- 4

  if(logscale) {
    if(percent==TRUE && theaxis=="x") {
      thelabels <- fmtRatioToPercent(axisrange) 
    }
    else if(percent==TRUE && theaxis=="y") {
      thelabels <- fmtRatioToPercent(axisrange)
    } 
    else {
      thelabels <- trimWhiteSpace(chopZeroes(fround(labelrange, digits)))
    }

    if(!is.null(zeroscore)) thelabels[1] <- "0"
  }
  else {
    if(!is.null(offset) || !is.null(zeroscore)) {
      stop(cgMessage("The offset and zeroscore arguments",
                     "must be set to NULL when",
                     "logscale=FALSE, as they are not compatible."))
    }
    if(difference==TRUE && (theaxis=="x" || theaxis=="y")) {
      thelabels <- fmtDifference(axisrange, digits=digits)
    }
    else {
      thelabels <- trimWhiteSpace(chopZeroes(fround(labelrange, digits)))
    }
  }

  if(maxonly) {
    thelabels <- c("", thelabels[2])
  }

  if(!grid) {
    if(theaxis=="y") {
      text(x=rep(par("usr")[1], 2), y=if(logscale) log10(axisrange) else axisrange,
           labels=paste("", thelabels), col="blue", adj=0, cex=cex)
    }
    else if(theaxis=="x") {
      text(y=rep(par("usr")[3], 2), x=if(logscale) log10(axisrange) else axisrange,
           labels=paste(thelabels,"\n"), col="blue", adj=0.5, cex=cex, xpd=NA)
    }
  }

  else { ## grid==TRUE
    if(theaxis=="y") {
      panel.text(x=rep(current.panel.limits()$xlim[1], 2),
                 y=if(logscale) log10(axisrange) else axisrange,
                 labels=paste("", thelabels), col="blue", adj=0, cex=cex)
    }
    else if(theaxis=="x") {
      panel.text(y=rep(current.panel.limits()$ylim[1], 2),
                 x=if(logscale) log10(axisrange) else axisrange,
                 labels=paste(thelabels, "\n"), col="blue", adj=0.5, cex=cex)
    }
  }

  invisible()
}

#####################

plotGrpNameTicks <- function(grpnames, settings, grid=FALSE) {
  ##
  ## PURPOSE:  Add group labels to x-axis
  ##
  padEnd <- function(x, length=2) {
    paste(x, paste(rep(" ", length), collapse=""), sep="")
  }
  
  numberofgrps <- length(grpnames)
  if(!grid) {
    if(settings$srt==0) {
      axis(1, at=1:numberofgrps, labels=grpnames,
           cex=settings$cex, adj=settings$adj)
    }
    else if(settings$srt!=0) {
      axis(1, at=1:numberofgrps, labels=FALSE)
      text(1:numberofgrps, par("usr")[3], srt = settings$srt,
           adj = settings$adj, cex=settings$cex,
           labels = padEnd(grpnames), ## to prevent axis overlap
           xpd = TRUE)
    }
  }
  else { ## grid, lattice
     panel.axis(side="bottom",
                at=1:length(grpnames),
                labels=if(settings$srt!=0) {
                  padEnd(grpnames)
                }
                else {
                  grpnames
                },
                tck=0, text.cex=settings$cex,
                rot=settings$srt,
                outside=TRUE)
  }
  invisible()

}


#####################

boxplotStamp <- function() {
  ##
  ## PURPOSE: Add explanatory message to boxplot rendering
  ##
  u <- par("usr")
  
  if(FALSE) {
    text(u[2] - 0.15 *(u[2] - u[1]),
       u[4] + 0.15 * (u[4] - u[3]),
       labels=paste("Open circles \"o\" ",
         "inside boxes are medians;",
         "Plus signs \"+\" are means", sep="\n"),
       cex=0.6, adj=0, xpd=TRUE)
  }

  mtext(text=paste("Open circles \"o\" ",
         "inside boxes are medians;",
         "Plus signs \"+\" are means", sep=" "),
        side=3,
        cex=0.6, col="gray", line=0.1, adj=0.5,
        xpd=NA)
  invisible()
}

#####################

errorBarGraphStamp <- function(mcadjust=FALSE, alphapercent, grid=FALSE) {
  ##
  ## PURPOSE: Add explanatory message to errorBarGraph rendering
  ##
  firstline <-
    if(!grid) {
      paste("Non-overlapping error bars indicate statistically\n",
            "significant difference(s)",
            sep=" ")
    }
    else {
      paste("Non-overlapping error bars within a panel indicate statistically\n",
            "significant difference(s) ",
            sep=" ")
    }
  
  secondline <- paste("at", alphapercent,
                      "% significance level",
                      if(mcadjust) {
                        "(Multiplicity adjusted)"
                      } else "", sep=" ")
  themessage <- paste(firstline,
                      secondline,
                      collapse="")

  if(!grid) {
    mtext(text=themessage, side=3, line=0.75, col="red", cex=0.7,
          xpd=NA, outer=FALSE, adj=0.5, font=2)
  }
  else {
    seekViewport(trellis.vpname("main"))
    grid.text(themessage,
              y=unit(-1, "lines"),
              gp=gpar(col="red", cex=0.7, fontface="bold"),
              just=c("centre","centre"))
  }
  
  invisible()
}

#####################

comparisonsGraphStamp <- function(mcadjust=FALSE, alphapercent, grid=FALSE, desc="") {
  ##
  ## PURPOSE: Add explanatory message to comparisonsGraph rendering
  ##
  themessage <- paste(
                      paste(" Error bars that do not cross the zero line indicate statistically\n",
                            "significant difference(s) at ",
                            sep=" "),
                      alphapercent,
                      "% significance level",
                      if(mcadjust) {
                        " (Multiplicity adjusted)"
                      } else "",
                      sep="")

  if(!grid) {
    mtext(text=themessage, side=3, line=0.75, col="red", cex=0.7,
          xpd=NA, outer=FALSE, adj=0, font=2)
    if(desc!="") {
      mtext(text=desc, side=3,
            cex=0.6, col="gray", line=0.1, adj=0.5,
            xpd=NA)
    }
    
  }
  else {
    seekViewport(trellis.vpname("main"))
    grid.text(themessage,
              y=unit(-1, "lines"),
              gp=gpar(col="red", cex=0.7, fontface="bold"),
              just=c("centre","centre"))
    if(desc!="") {
      grid.text(desc,
                x=unit(0.5, "npc"),
                y=unit(-2.75, "lines"),
                gp=gpar(col="gray", cex=0.6))
    }
  }
  
  invisible()
}


#####################

errorBarGraphApproximateStamp <- function(grid=FALSE, msg=NULL) {
  ##
  ## PURPOSE: Add warning message about unequal sample sizes
  ##  and/or standard error estimates
  ##
  if(is.null(msg)) msg <- "Error bars are approximate."
  if(!grid) {
    if(msg!="") {
      mtext(text=msg, side=3,
            cex=0.6, col="gray", line=0.1, adj=0.5,
            xpd=NA)
    }
  }
  else {
    seekViewport(trellis.vpname("main"))
    if(msg!="") {
      grid.text(msg,
                x=unit(0.5, "npc"),
                y=unit(-2.75, "lines"),
                gp=gpar(col="gray", cex=0.6))
    }
  }

  invisible()
}

#####################

trimWhiteSpace <- function(x) {
  ##
  ## PURPOSE:   Laziness motivated wrapper to
  ## trim leading and trailing spaces in a character string
  ## NOTE: x can be a vector of strings, and empty contents
  ## will be converted to an empty string.
  ## EXAMPLE USAGE: trimWhiteSpace(c(" 37 ","abc ","jfk", 4576))
  ##
  if(length(x)==0) return("")
  y <- gsub("^\\s+","", x, perl=TRUE) ## leading
  return(gsub("\\s+$","", y, perl=TRUE)) ## trailing
}

#####################

chopZeroes <- function(x) {
  ##
  ## PURPOSE: For tickmark labels, delete any superfluous trailing zeroes after a
  ## decimal point.  In effect, convert the character item to look like an integer.
  ## EXAMPLE USAGE: chopZeroes(c("30","12.0","33.000","478.020","65.20","0.1205"))
  ##
  if(length(x)==0) return("")
  gsub("\\.0*0+$", "", as.character(x), perl=TRUE)
}

#####################

fmtRatioToPercent <- function(x, digitsbelow10=1) {
  ##
  ## PURPOSE: Map ratios to percent units for tick labels
  ##
  y <- 100 * (x - 1)
  outfmt <- trimWhiteSpace(ifelse(abs(y) < 9.5 & y!=0, formatC(round(y, digitsbelow10),
                                                   format="f", digits=digitsbelow10),
                                  chopZeroes(round(y, 0))))
  outfmt[as.numeric(outfmt)==0] <- "0"
  return(outfmt)
}

#####################

fmtDifference <- function(x, digits) {
  ##
  ## PURPOSE: Prettify Differences for tick labels 
  ##
  trimWhiteSpace(ifelse(abs(x) < 10,
                        chopZeroes(formatC(round(x, digits=digits),
                                           format="f", digits=digits)),
                        chopZeroes(round(x, 0))))
}

#####################

fmtRatio <- function(x) {
  ##
  ## PURPOSE: Prettify Ratios for tick labels 
  ##
  trimWhiteSpace(ifelse(abs(x) < 0.1, formatC(round(x, 3), format="f", digits=3),
                        chopZeroes(round(x, 2))))
}

#####################

fmtPercent <- function(x, decimalthreshold=10) {
  out <- fround(x, 0)
  for(i in seq(along=out)) {
    if(out[i]=="-100") {
      out[i] <- "< -99.9" 
    }
    ## else if(out[i]=="0") {
    else if(abs(x[i]) < decimalthreshold) {
      out[i] <- fround(x[i], 1)
    }
  }
  return(out)
}

#####################

fmtPvalue <- function(x) {
  out <- fround(x, 3)
  for(i in seq(along=x)) {
    if(x[i] < 0.0005) {
      out[i] <- "<0.001" 
    }
  }
  return(out)
}


#####################

cgDevice <- function(cgtheme=TRUE, new=FALSE, ...) {
  ##
  ## PURPOSE: Ensure trellis device with limited color scheme
  if(cgtheme) {
    thetheme <- list(background = list(col = "white"),
                     strip.shingle = list(col = "white"),
                     strip.background = list(col = "white"))
  }
  else {   
    thetheme <- trellis.par.get()
  }
  if(.Device=="null device") {
    new <- TRUE
  }
  trellis.device(theme=thetheme, new=new, ...)
  invisible()
}

#####################

contrastMatrix <- function(listoflevels,
                           type = c("allgroupstocontrol",
                             "pairwisereflect",
                             "pairwise",
                             "successivetimes")
                           ) {
  ##
  ## PURPOSE: Utility function to construct comparisons specifications
  ## NOTE: The contrast matrix generation code is borrowed and adapted from the
  ## contrMat function of the R 'multcomp' package
  ##
  if(mode(listoflevels) != "list") {
    listoflevels <- list(listoflevels)
  }

  ## start of R multcomp::contrMat code adaptation
  type <- match.arg(type) ## ensure default is 'allgroupstocontrol'
  
  ## Note that "grpnames" is just a generic term for
  ## the different levels of a factor
  genAtOneLevel <- function(grpnames, type=type) {
    CM <- c()
    rnames <- c()
    k <- length(grpnames)
    kindx <- 1:k
    
    switch(type,
           allgroupstocontrol = {
             for (i in 2:k) CM <- rbind(CM, as.numeric(kindx == i) - 
                                        as.numeric(kindx == 1))
             rnames <- paste(grpnames[2:k], " vs. ", grpnames[1], sep = "")
           },
           pairwisereflect = {
             for (i in 1:k) {
               for (j in (1:k)[-i]) {
                 CM <- rbind(CM, as.numeric(kindx == j) - as.numeric(kindx == 
                                              i))
                 rnames <- c(rnames, paste(grpnames[j], " vs. ", grpnames[i], 
                                           sep = ""))
               }
             }
           },
           pairwise = {
             for (i in 1:(k - 1)) {
               for (j in (i + 1):k) {
                 CM <- rbind(CM, as.numeric(kindx == j) - as.numeric(kindx == 
                                              i))
                 rnames <- c(rnames, paste(grpnames[j], " vs. ", grpnames[i], 
                                           sep = ""))
               }
             }
           },
           successivetimes = {
             for(i in 2:k) {
               CM <- rbind(CM, as.numeric(kindx == i) - 
                           as.numeric(kindx == (i-1)))
               rnames <- c(rnames, paste(grpnames[i], " vs. ", grpnames[i-1], sep = ""))
             }
           }
           )
    dimnames(CM) <- list(rnames, grpnames)
    CM

  }
  
  ## end of R multcomp::contrMat code adaptation
  listLs <- lapply(listoflevels, genAtOneLevel, type=type)
  if(length(listLs)==1) {
    Lcomparisonmat <- listLs[[1]]
  }
  else {
    Lcomparisonmat <- blockDiag(listLs)
  }
  
  return(Lcomparisonmat)
  
}

#####################

blockDiag <- function(listmat) {
  ##
  ## PURPOSE: Create a block diagonal matrix from list of matrices
  ## NOTE: adapted from Ben Bolker contribution on R-help
  ## http://www.r-project.org/nocvs/mail/r-help/2002/1533.html
  ##
  nr <- sapply(listmat,nrow)
  nc <- sapply(listmat,ncol)
  cumnc <- cumsum(nc)
  NR <- sum(nr)
  NC <- sum(nc)
  rowfun <- function(m,zbefore,zafter) {
    cbind(matrix(0,ncol=zbefore,nrow=nrow(m)),m,
          matrix(0,ncol=zafter,nrow=nrow(m)))
  }
  ret <- rowfun(listmat[[1]],0,NC-ncol(listmat[[1]]))
  for (i in 2:length(listmat)) {
    ret <- rbind(ret,rowfun(listmat[[i]],cumnc[i-1],NC-cumnc[i]))
  }
  dimnames(ret)[[2]] <- as.vector(unlist(lapply(listmat, function(x)
                                                dimnames(x)[[2]])))
  attr(ret, "colsperblock") <- nc
  attr(ret, "rowsperblock") <- nr
  ret
} 

#####################

rangeExtend <- function(x, pct=4) {
  ##
  ## PURPOSE: Utility function to extend axis ranges, as sometimes trellis
  ## or lattice/grid functions in R seem to come a bit short, and points or
  ## labels spill over the plot region boundaries
  ##
  if(length(pct)==1) {
    pct <- list(minside=4, maxside=4)
  }
  rangelength <- diff(range(x, na.rm=TRUE))
  minx <- max(x, na.rm=TRUE) - (1+(pct$minside/100))*rangelength
  maxx <- min(x, na.rm=TRUE) + (1+(pct$maxside/100))*rangelength
  c(minx, maxx)
}

#####################

getNumDigits <- function(dfr, max=4) {
  ##
  ## PURPOSE: Utility function to determine the maximum number of digits
  ## used in values of an input set of data
  ## argument dfr can be a data frame, list, matrix, or atomic vector
  ##
  thedata <- as.character(as.vector(unlist(dfr)))
  maxdigits <- max(nchar(gsub("(\\d*)(\\.*)(\\d*)", "\\3", thedata,
                              perl=TRUE)))
  if(maxdigits > 4) maxdigits <- max
  return(maxdigits)
}

#####################

makeCensored <- function(x, boundary, direction=">") {
  fmtx <- x
  for(i in seq(along=x)) {
    if(x[i] >= boundary) {
      fmtx[i] <- paste(direction, x[i], sep="")
    }
  }
  return(fmtx)
}


#####################

cgLineColors <- c("black", "blue", "green", "red", "orange", "brown", "yellow",
                  "darkblue", "darkgreen", "darkgrey",
                  "black", "blue", "green", "red", "orange", "darkgrey") 

#####################


multcompInform <- function() {
  cat(cgMessage("Some time may be needed as the",
                "critical point from the multcomp::summary.glht",
                "function call is calculated. Please wait...\n\n",
                warning=TRUE))
  flush.console()
  invisible()
}

multcompDone <- function(model) {
  cat(cgMessage(paste("...Done. Critical point from", model,"fit",
                      "is calculated.\n\n"),
                warning=TRUE))
  flush.console()
  invisible()
}

#####################

isAllEqual <- function(x) {
  ## x is assumed to be an atomic vector
  return(length(unique(x))==1)
}

#####################

makeEndptLabel <- function(endptname, endptunits) {
  ## set up for input into expression
  e1c <- if(is.expression(endptname)) {
    paste(endptname) }
  else {
    paste("plain('",endptname,"')", sep="")
  }
  e2c <-
    if(is.expression(endptunits)) {
      paste(endptunits) }
    else {
      paste("plain('",endptunits,"')", sep="")
    }
  noe2c <- (e2c=="plain('')")
  expr <- parse(text=paste(
                  "paste(",
                  e1c,
                  if(!noe2c) " , ",
                  if(!noe2c) "' (',",
                  if(!noe2c) e2c,
                  if(!noe2c) " , ",
                  if(!noe2c) "')'",
                  ")",
                  sep=" "))
  return(expr)
}

#####################

catCharExpr <- function(char, expr, rev=FALSE) {
  if(!rev) {
    out <- parse(text=paste("paste(",
                   paste("plain('",char,"')", sep=""),
                   " , ",
                   paste(expr),
                   ")",
                   sep=" "))
  }
  else {
    out <- parse(text=paste("paste(",
                   paste(expr),
                   " , ",
                   paste("plain('",char,"')", sep=""),
                   ")",
                   sep=" ")) 

  }
  return(out)
}



#####################

residualgrptrend.helper <- function(tukeyresids, fitteds, tukeyresids2=NULL,
                                    status=NULL, grpf, desc="") {
  ## 
  ## PURPOSE: Determine a trend for residuals (y) against
  ## fitted values (x), with likely option of rescaling the fitted values
  ## to integer values afterwards: 1 to number of groups.
  ##
  trend.index <- order(fitteds)
  
  thegrpf <- factorInSeq(grpf[trend.index])
  tukeyresids <- tukeyresids[trend.index]
  fitteds <- fitteds[trend.index]

  if(!is.null(status)) {
    ## Assume proper censored data format (survival:::Surv())
    ## no checking
    has.censored <- TRUE
    tukeyresids2 <- tukeyresids2[trend.index]
    status <- status[trend.index]
    thesurvobject <- survival::Surv(time=tukeyresids, time2=tukeyresids2,
                                    event=status, type="interval")
    if(nlevels(thegrpf) < 7) {
      ## Fit "group means" nominal AFT-normal model
      thefit <- try(survreg(thesurvobject ~ -1 + thegrpf,
                            dist="gaussian", maxiter=100))
    }
    else {
      ## Need different data format for VGAM:::vgam()
      forvgam.dfr <- data.frame(x=fitteds, y=tukeyresids)
      extra <- list(leftcensored = (status==2),
                    rightcensored = (status==0))
      thefit <- try(VGAM::vgam(y ~ s(x, 4), VGAM::cennormal(zero=2),
                               data=forvgam.dfr, trace=FALSE,
                               extra=extra))
      if(inherits(thefit, "try-error")) {
        thefit <- try(VGAM::vgam(y ~ s(x, 3), VGAM::cennormal(zero=2),
                                  data=forvgam.dfr, trace=FALSE,
                                  extra=extra))
      }
      if(inherits(thefit, "try-error")) {
        thefit <- try(VGAM::vgam(y ~ s(x, 2), VGAM::cennormal(zero=2),
                                  data=forvgam.dfr, trace=FALSE,
                                  extra=extra))
      }
      else if(inherits(thefit, "try-error")) {
        ## Fit "group means" nominal AFT-normal model if smoothing fails
        thefit <- try(survreg(thesurvobject ~ -1 + fitteds,
                              dist="gaussian", maxiter=100))
      }
    }
    trendtukeyresids <- if(class(thefit)[1]=="survreg") {
      predict(thefit) } else if(class(thefit)[1]=="vgam") { predict(thefit)[, 1] }
  }
  else {
    has.censored <- FALSE
    trendtukeyresids <- lowess(fitteds, tukeyresids)$y
  }
  
  dfr <- data.frame(grpc=I(as.character(thegrpf)),
                    grpn=as.numeric(thegrpf),
                    tukeyresids=tukeyresids,
                    trendtukeyresids=trendtukeyresids,
                    type=rep(desc, length(trend.index)))
  if(has.censored) {
    dfr$status <- status
  }
  return(dfr)
}

#####################


fround <- function(x, digits=0) {
  format(round(x, digits=digits), nsmall=digits)
}

######################

fround.charcens <- function(x, digits=0) {
  ## generalization when x has "<" and ">" representations
  for(i in seq(along=x)) {
    if(is.na(x[i]) || x[i]=="<NA>") next  
    else if(regexpr("<|>", x[i]) > 0) {
      pieces <- strsplit(x[i], "<|>")
      x[i] <- paste(substr(x[i], 1, 1), fround(as.numeric(pieces[[1]][2]), digits),
                    sep="")
    }
    else {
      x[i] <- fround(as.numeric(x[i]), digits)
    }
  }
  return(x)
}


#####################

chop.matrix <- function(x) {
  n <- nrow(x)
  return(x[1:(n-1), 1:(n-1)])
}


#####################

prepare <- function(type, ... ) {
  ##
  ## PURPOSE: Read-in the input data frame, and create
  ## an object with it 
  ## and add metadata for cg methods.
  ## Essentially supplies alias to specific prepare* methods
  type <- validArgMatch(type, c("onefactor","unpairedgroups",
                                "pairedgroups","paireddifference",
                                "repmes","repeatedmeasures"))

  if(is.element(type, c("onefactor","unpairedgroups"))) {
    return(switch(type,
                  "onefactor"={ prepareCGOneFactorData(...) },
                  "unpairedgroups"={ prepareCGOneFactorData(...) }
    ))
  }
  else if(is.element(type, c("paireddifference","pairedgroups"))) {
    return(switch(type,
                  "paireddifference"={ prepareCGPairedDifferenceData(...) },
                  "pairedgroups"={ prepareCGPairedDifferenceData(...) }
    ))
  }
   
  invisible()

}

#####################


stripmiss <- function(x) x[!is.na(x)]

#####################

makeTickMarks <- function(ticklabels, tickmarks, percent=FALSE,
                          offset=NULL) {
  validList(ticklabels,
            names=c("mod","marks"),
            argname="ticklabels")
  validBoolean(percent)
  mod <- eval(parsePartialName(names(ticklabels), "mod",
                               prefix="ticklabels$"))
  mod <- validArgMatch(mod, choices=c("add", "replace"))
  marks <- eval(parsePartialName(names(ticklabels), "marks",
                                 prefix="ticklabels$"))
  
  if(mod=="replace") {
    tickmarks <- as.numeric(marks)
    if(percent) { tickmarks <- pctToRatio(tickmarks) }
    if(!is.null(offset)) tickmarks <- (tickmarks + offset)
    names(tickmarks) <- marks
  }
  else if(mod=="add") {
    addedtickmarks <- as.numeric(marks)
    if(percent) { addedtickmarks <- pctToRatio(addedtickmarks) }
    if(!is.null(offset)) addedtickmarks <- (addedtickmarks + offset)
    names(addedtickmarks) <- marks
    tickmarks <- sort(c(tickmarks, addedtickmarks))
  }
  
  return(tickmarks)
}

#####################

scaleVar <- function(x, endptscale, logbase=10, percent=FALSE) {
  ##
  ## PURPOSE: Simple function to help work with variables
  ## that may or may not be transformed
  endptscale <- validArgMatch(endptscale, c("log","original"))
  validBoolean(percent)
  if(endptscale=="log") {
    if(percent) { x <- pctToRatio(x) }
    return(log(x, logbase))
  }
  else return(x)
}

#####################

makeContrastVec <- function(posplaces=NULL, negplaces=NULL, length,
                            posvalues=rep(1, length(posplaces)),
                            negvalues=rep(-1, length(negplaces))
                            ) {
  ##
  ## PURPOSE: Simple function to create vectors in
  ## a contrast matrix for custom comparisons

  ## WARNING: Argument checking needs to be implemented!
  vec <- numeric(length) ## initialization
  vec[posplaces] <- posvalues
  vec[negplaces] <- negvalues

  return(vec)
}

#####################

orderPairedGroups <- function(x, refgrp) {
  ## PURPOSE: For the PairedDifference
  ## methods, order the 2 group names
  ## such that the refgrp is first and
  ## derive the difference between then
  if(x[1]==refgrp) {
    return(x)
  }
  else {
    c(x[2], x[1])
  }
}



###########################

qminmin <- function(y, x, q) {
  ## PURPOSE: compute quantiles from survival curve estimate
  ##          to aid in boxplot graphing
  ## NOTE: adapted from minmin helper function within
  ## survmean helper function from the survival package.
  ## (The minmin function used to be in the print.survfit function
  ## in earlier survival package versions).
  p <- 1 - q
  keep <- (!is.na(y) & y <= p)
  if (!any(keep)) 
    NA
  else {
    x <- x[keep]
    y <- y[keep]
    if (y[1] == p && any(y < p)) 
      (x[1] + x[min(which(y < p))])/2
    else x[1]
  }
}

######################
