##*****************************************************************************
##      DO NOT ROXYGENISE THIS!!!
## 
##' ##' Temporal aggregation of a Marked Process, leading to block maxima
##' or \eqn{r}-largest observations.
##'
##' The data frame given in \code{OTdata} contains the \emph{events}
##' (or \emph{arrivals}) given by the \code{date} column, as well as
##' one mark column. Depending on the argument \code{MAX.r}, the
##' maxima or the \eqn{r}-largest observations of the marks is
##' computed for each time block. When known gaps exist in the data
##' and when they are given in \code{OTmissing}, a block for which the
##' total duration of gaps is too large will be omitted.
##'
##' @title Temporal aggregation of a Marked Process
##'
##' @param OTdata Data frame containing a \code{POSIXct} column
##' \code{date} and the marks variable.
##' 
##' @param OTmissing Optional data frame with columns \code{start} and
##' \code{end} (coerced to \code{POSIXct}) giving the beginning and the
##' end of gaps.
##' 
##' @param start An object coerced to POSIXct indicating the begining
##' of reliable/usable information. Unless this is a begining of block
##' (1-st of january for years), the 1-st block will begin \emph{after}
##' \code{start} in order to use only qualified information.
##' 
##' @param end An object indicating the end of the reliable/usable
##' information. Unless this is a end of block
##' (1-st of january for years), the last block will end \emph{before}
##' \code{end} in order to use only qualified information. 
##'
##' @param MAX.r \emph{Target} number of observations in the
##' blocks. Can be of length 1 (same number of observations for all
##' blocks) or of length equal to the number of blocks, the values
##' being then for the blocks in the same order. In both cases, the
##' target number may be impossible to reach because of a smaller
##' number of events in the block. If \code{infMAX} is \code{TRUE} the
##' target number of observations will be reached by filling if needed
##' with \code{-Inf} values. The rationale for this is that a
##' non-existing event is assumed to have an arbitrarily small
##' mark.
##' 
##' @param blockDuration Duration of the blocks. Can only be
##' \code{"year"} for now.
##'
##' @param monthGapStat Logical. Setting it to \code{TRUE} will
##' compute statistics concerning the gaps and return them or show them
##' in a plot.
##' 
##' @param maxMissingFrac Maximal fraction of a block duration
##' (between 0 and 1) that can be missing without leading to a
##' \code{NA} aggregated value.
##'
##' @param dataFrames If \code{TRUE}, the result will contain data
##' frames as found similar to those found in an object with class
##' \code{"Rendata"}. If \code{FALSE} the result will contain
##' \emph{list} and \emph{vector} objects, similar to those used as
##' inputs in the \code{\link{Renouv}} function under the names
##' \code{MAX.data} and \code{MAX.effDuration}. Note however, that
##' \code{-Inf} values can be found in these objects when
##' \code{infMAX} is \code{TRUE}.
##'
##' @param infMAX If \code{FALSE}, the target number of values the
##' blocks will generally not be reached, because the total number of
##' events in a block can be lower than the target number. Then, the
##' target number value is revised to the number of found values in
##' each block. If \code{TRUE}, the target number of values is reached
##' by filling the values with \code{-Inf} and the datetimes with
##' (\code{POSIXct}) \code{NA}s.
##'
##' @param plot If \code{TRUE} a simple plot is shown.
##'
##' @param plotType Character controling the plot. With
##' \code{"max"}, the block maxima are shown. With \code{"gap"}
##' the daily and monthly gap rates are shown. This is possible
##' when suitable information concerning gaps is provided in
##' \code{OTmissing}.
##' 
##' @param jitterSeed Random seed for jittering. Used only when
##' \code{plot} is \code{TRUE}, \code{plotType} is \code{"gap"}
##' and when suitable information is provided in \code{OTmissing}.
##' 
##' @param trace Integer level of verbosity.
##'
##' @param ... Other arguments to be passed to \code{plot}.
##'
##' @return A list, the content of which depends on the value of
##' \code{dataFrames}. If this value is \code{TRUE}, the following
##' elements are returned.
##'
##' \item{MAXdata}{
##' 
##' A data frame of largest values by block with one row for each
##' observation. The largest values are given as columns with names
##' equal to those in the \code{OTdata} data frame.
##'
##' }
##'
##' \item{MAXinfo}{
##'
##' A data frame describing the blocks, with one row by block. The two
##' (\code{POSIXct}) columns \code{"start"} and \code{"end"} provide
##' the beginning and the end of the block.  The numeric column
##' \code{duration} gives the \emph{effective duration} (in year)
##' within block.
##'
##' }
##'
##' \item{probMissing}{
##'
##' A vector with values corresponding to the days in a block
##' (year). Each value is a estimation of the probability that the day
##' falls in a gap.
##'
##' }
##'
##' If \code{dataFrames} is \code{FALSE}, the list still contains
##' \code{probMissing} as before, as well as other lists as used in
##' \code{\link{Renouv}}.
##' 
##'
##' \item{effDuration, r}{
##'
##'  Vectors containing the effective duration (\emph{in years}) and
##' number of value for the blocks.
##'  
##' }
##' 
##' \item{data}{
##'
##'  List of maxima or \eqn{r}-largest values for the blocks.
##'  
##' }
##'
##' \item{monthGapStat, monthGapTS}{
##'
##' Summary information concerning
##' gaps, \code{monthGapStat} is \code{TRUE} and if relevant
##' information is provide via the the \code{OTmissing} formal.
##' 
##' }
##' 
##' @note Remind that even when \code{maxMissingFrac} is set to its
##' maximum value 1.0, there can still be blocks with no data. When
##' the result is intended to be used in the \code{\link{Renouv}}
##' function, the formal \code{dataFrames} should be \code{FALSE}; the
##' elements \code{data} and \code{effDuration} can then be passed as
##' \code{MAX.data} and \code{MAX.effDuration}. At the time
##' \code{infMAX} should also then be set to \code{FALSE} since
##' \code{-Inf} values are not yet allowed in the \eqn{r}-largest
##' values.
##'
##' @caution This function is intended to be used with a time
##' precision limited to about one day.
##'
##' @examples
##' ## use Dunkerque data
##' OTdata <- Dunkerque$OTdata; OTmissing <- Dunkerque$OTmissing
##' ## allow up to 50\% gap in a block, or only 5\%
##' MAX1 <- OT2MAX(OTdata = OTdata, OTmissing = OTmissing,
##'                maxMissingFrac = 0.5,
##'                main = "impact of the 'maxMissingFrac' formal")
##' MAX2 <- OT2MAX(OTdata = OTdata, OTmissing =OTmissing,
##'                prefix = "Max", maxMissingFrac = 0.05, plot = FALSE) 
##' lines(MAX2$MAXdata$date, MAX2$MAXdata$Surge, type = "h", col = "red", lwd = 3)
##' legend("topleft", lw = c(1, 3), col = c("black", "orangered"),
##'        legend = c("50\% max", " 5\% max"))
##' 
##' ## plot the probability of being in a gap along a year
##' plot(MAX2$probMissing, type ="l", ylim = c(0, 1),
##'      xlab = "day in year", ylab = "prob")
##'
##' ## r-largest obs for r = 4
##' MAX3 <- OT2MAX(OTdata, OTmissing = OTmissing, MAX.r = 4,
##'                maxMissingFrac = 0.9, infMAX = TRUE,
##'                dataFrames = FALSE, trace = TRUE)
##'
##' ## plot summary information concerning gaps.
##' MAX4 <- OT2MAX(OTdata = OTdata, OTmissing = OTmissing,
##'                maxMissingFrac = 0.5,
##'                main = "impact of the 'maxMissingFrac' formal",
##'                plotType = "gap")
##' 
##' ## time series plot (only <= 10 months)
##' plot(MAX4$monthGapTS[ , c(1:4)])
##' 
##' ## much better 
##' \dontrun{
##'     require(lattice)
##'     xyplot(MAX4$monthGapTS)
##' }
##' 
OT2MAX <- function(OTdata,
                   OTmissing = NULL,
                   start = NULL,
                   end = NULL,
                   MAX.r = 1L,
                   blockDuration = "year",
                   monthGapStat = TRUE,
                   ## breaks = NULL, blockNames = NULL,
                   maxMissingFrac = 0.05,
                   dataFrames = FALSE,
                   infMAX = FALSE,
                   plot = TRUE,
                   plotType = c("max", "gaps"),
                   jitterSeed = 123,
                   trace = 0L,
                   ...) {
    
    plotType <- match.arg(plotType)
    
    ## the time zone is given by 'OTdata$date'
    if (is(OTdata$date, "POSIXct")) {
        TZ <- attr(OTdata$date, "tz")
    } else {
        warning("'OTdata$date' is not of class \"POSIXct\". Coerced ",
                " with time zone \"GMT\"")
        TZ <- "GMT"
        OTdata$date <- as.POSIXct(OTdata$date, tz = "TZ")
    }

    ## build the sequence of 366 days
    unit <- blockDuration
    Year <- seq(from = as.POSIXct("2000-01-01", tz = TZ),
                to = as.POSIXct("2000-12-31", tz = TZ),
                by = "day")
    labYear <- format(Year, "%m-%d")
    
    if ( !is.data.frame(OTdata) || ! ("date" %in% colnames(OTdata)) ) {
        stop("'OTdata' must be a data frame object with a 'date' column")
    }

    ## 
    nData <- nrow(OTdata)
    
    ## find variable name
    vars <- names(OTdata)
    vars <- vars[is.na(match(vars, c("date", "comment")))]
    if (length(var) > 1L) stop("variable name could not be guessed")
    
    ## check 'OTmissing', if provided
    if (!is.null(OTmissing)) {
        if ( !is.data.frame(OTmissing) || !("start" %in% colnames(OTmissing)) ||
            !("end" %in% colnames(OTmissing)) ) {
            stop("'OTmissing' must be a data frame object 'start' ",
                 "and 'end' columns")
        }
        miss <- TRUE
        nMP <- nrow(OTmissing)
        if (!is(OTmissing$start, "POSIXct")){
            OTmissing$start <- as.POSIXct(OTmissing$start, tz = TZ)
        } else {
            if (attr(OTmissing$start, "tz") != TZ) {
                warning("the timezone of 'OTmissing$start' does not match ",
                        "that of 'OTdata$date")
            }
        }
        if (!is(OTmissing$end, "POSIXct")){
            OTmissing$end <- as.POSIXct(OTmissing$end, tz = TZ)
        } else {
            if (attr(OTmissing$end, "tz") != TZ) {
                warning("the timezone of 'OTmissing$end' does not match ",
                        "that of 'OTdata$date")
            }
        }
    } else {
        miss <- FALSE
        if (monthGapStat) {
            monthGapStat <- FALSE
            warning("'monthGapStat' forced to FALSE since no information was ",
                    "provided in 'OTmissing'")
        }
    }
    
    ## compute the begining/end of the data
    if (is.null(start)) {
        if (miss) {
            if (OTmissing$start[1L] < OTdata$date[1L]) {
                start <- OTmissing$start[1L]
            } else start <- OTdata$date[1L]
        } else start <- OTdata$date[1L]
    } else {
        if (!is(start, "POSIXct")) start <- as.POSIXct(start, tz = TZ)
        else if (attr(start) != TZ) {
            warning("the timezone of 'start' does not match ",
                    "that of 'OTdata$date")
        }
    }
    if (is.null(end)) {
        if (miss) {
            if (OTmissing$end[nMP] > OTdata$date[nData]){
                end <- OTmissing$end[nMP]
            } else end <- OTdata$date[nData]
        } else end <- OTdata$date[nData]
    } else {
        if (!is(end, "POSIXct")) end <- as.POSIXct(end, tz = TZ)
        else if (attr(end) != TZ) {
            warning("the timezone of 'end' does not match ",
                        "that of 'OTdata$date")
            
        }
    }
    
    ## compute the begining/end of the blocks. This is obtained by SHRINKING
    ## the period.
    startBreaks <- start
    endBreaks <- end
    startBreaksTest <- as.POSIXct(format(startBreaks, "%Y-01-01"), tz = TZ)
    if (startBreaks != startBreaksTest) {
        year <- as.integer(format(startBreaks, "%Y")) + 1
        startBreaks <- as.POSIXct(paste(year, "-01-01 00:00", sep = ""), tz = TZ)
    }
    endBreaks <- as.POSIXct(format(endBreaks, "%Y-01-01"), tz = TZ)
    ## if (endBreaks != endBreaksTest) {
    ##    year <- as.integer(format(endBreaks, "%Y")) - 1
    ##    endBreaks <- as.POSIXct(paste(year, "-01-01 00:00", sep = ""), tz = TZ)
    ##}
    
    ## find blocks and interpolation grid for 'probMissing'
    Breaks <- seq(from = startBreaks, to = endBreaks, by = unit)
    nBreaks <- length(Breaks)
    
    blockNames <- format(Breaks[-nBreaks], "%Y")
    Block <- cut(x = OTdata$date, breaks = Breaks)
    levels(Block) <- blockNames
    nBlock <- nlevels(Block)
    blockDur <- as.numeric(diff(Breaks), units = "days")
    Breaks2 <- seq(from = startBreaks, to =  endBreaks, by = "day")
    
    if (length(MAX.r) != 1L && length(MAX.r) != nBlock) {
        stop("the length of 'MAX.r' must be 1 or the number of blocks") 
    }
    MAX.r <- rep(MAX.r, length.out = nBlock)
    names(MAX.r) <- blockNames
    if (trace) {
        cat("Number of events by block \n")
        print(tapply(OTdata[ , vars], Block, length))
    }
    
    ##===========================================================================
    ## Compute the cumulated duration of the gaps from the orgin and
    ## perform a linear intrpolation at 'Breaks' with result in
    ## 'interp'. The difference between the values at two successive
    ## breaks (i.e. the begining and the end of one block) gives the gap
    ## duration in the block. This is stored in the 'MAXinfo' data frame.
    ##
    ## The 'interp2' interpolation is done at a fine grid level (days) and
    ## is used to compute the probability of being in a gap along a year.
    ##
    ## CAUTION : remind that the concatenation method 'c' used with POSICXct
    ## objects loose the "tz" attribute"
    ##
    ##     c(as.POSIXct("2015-01-01", tz = "ACT"),
    ##       as.POSIXct("2015-01-20", tz = "ACT"))
    ##
    ## is NOT the same thing as
    ##
    ##     as.POSIXct(c("2015-01-01", "2015-01-20"), tz = "ACT")
    ##
    ## so we must use data frames and 'rbind'.
    ## ==========================================================================
    if (!miss) {
        ms <- NULL
        mTS <- NULL
        blockDurMiss <- rep(0.0, nBlock)
        keepBlock <- rep(TRUE, nBlock)
    } else {
        ## put 'start' and 'end' of missing periods in the right order
        ## gap = 1 means that a gap follows the date, val = 0 means
        ## that a non-gap follows.
        dfMiss <- data.frame(dt = OTmissing$start,
                             gap = rep(1, length(OTmissing$start)))
        dfMiss <- rbind(dfMiss,
                        data.frame(dt = OTmissing$end,
                                   gap = rep(0, length(OTmissing$end))))
        
        dfMiss <- dfMiss[order(dfMiss$dt), ]
        ## add the cumsum column
        cum <- cumsum(c(0.0, as.numeric(diff(dfMiss$dt), units = "days")) *
                          (1.0 - dfMiss$gap))
        dfMiss <- cbind(dfMiss, cum = cum)
        ## add a row before and/or after if needed.
        if (dfMiss$dt[1L] > Breaks[1L]) {
            dfMiss <- rbind(data.frame(dt = Breaks[1L], gap = 0, cum = 0),
                            dfMiss)
        }
        nMiss <- nrow(dfMiss)
        if (dfMiss$dt[nMiss] < Breaks[nBreaks]) {
            dfMiss <- rbind(dfMiss,
                            data.frame(dt = Breaks[nBreaks],
                                       gap = 0,
                                       cum = dfMiss$cum[nMiss]))
        }
        ## interpolate: linear interpolation for cumulative duration,
        ## constant for the state (in Gap or not)
        interp <- approx(x = dfMiss$dt, y = dfMiss$cum, xout = Breaks,
                         method = "linear")
        blockDurMiss <- as.numeric(diff(interp$y), units = "days")
        keepBlock  <- ( (blockDurMiss / blockDur) < maxMissingFrac )
        names(keepBlock) <- levels(Block)
        BlockMiss <- as.integer(cut(x = dfMiss$dt, breaks = Breaks))
        interp2 <- approx(x = dfMiss$dt, y = dfMiss$gap, xout = Breaks2,
                          method = "constant")
        ## if monthly stats are required, interpolate at month breaks
        if (monthGapStat) {
              monthBreaks <- seq(from = startBreaks,
                                 to = endBreaks,
                                 by = "month")
              monthDur <- round(as.numeric(diff(monthBreaks), units = "days"))
              monthInterp <- approx(x = dfMiss$dt, y = dfMiss$cum,
                                    xout = monthBreaks,
                                    method = "linear")
              monthBlockDurMiss <- as.numeric(diff(monthInterp$y), units = "days") 
              nM <- length(monthBreaks)
              monthStart <- seq(from = as.POSIXct("2001-01-01", tz = TZ),
                                by = "month",
                                length.out = 13)
              monthLab <- format(monthStart[1:12], "%b")
              month <- factor(format(monthBreaks[-nM], "%b"), levels = monthLab)
              
              ## monthBlockNames <- format(monthBreaks, "%Y-%B")
              ms <- data.frame(start = monthBreaks[-nM],
                               end = monthBreaks[-1],
                               year =  as.numeric(format(monthBreaks[-nM], "%Y")),
                               ## month = as.numeric(format(monthBreaks[-nM], "%m")),
                               month = month,
                               ## monthName = format(monthBreaks[-nM], "%B"),
                               dur = monthDur,
                               durMissing = monthBlockDurMiss,
                               missing = monthBlockDurMiss / monthDur)
              rownames(ms) <- format(monthBreaks[-nM], "%Y-%m")
              
              mTS <- tapply(ms$missing, list(ms$year, ms$month), mean)
              mTS <- ts(mTS, start = min(ms$year))
              
          } else {
              ms <- NULL
              mTS <- NULL
          }
    }
    
    effDur <- 1 - (blockDurMiss / blockDur)
    effDur[abs(effDur) < 0.001] <- 0.0
    MAXinfo <- data.frame(start =  Breaks[-nBreaks],
                          end = Breaks[-1L],
                          duration = effDur)
    rownames(MAXinfo) <- blockNames
    
    ##===========================================================================
    ## construct 'probMiss' and fill 'MAXdata'. The 'MAXdata' object is
    ## in a first time as a list, and then is coerced to data frame.
    ##
    ## CAUTION: when coetcing 'Block' to a factor, remind that some year
    ## may have no data and then will need care to appear in the levels. In
    ## order to get the right levels, 'block' is first considered as a character
    ## vector (and not as an integer vector).
    ## ==========================================================================
    probMiss <- rep(0, 365)    
    first <- TRUE
    
    for (i in 1L:nBlock) {
        
        ind <- (OTdata$date >= Breaks[i]) & (OTdata$date <= Breaks[i + 1L])
        
        if (keepBlock[i]) {
            ## if there is at least an event
            if (any(ind)) {
                ri <- sum(ind)
                Ti <- OTdata[ind, "date"]
                zi <- OTdata[ind, vars]
                indi <- rev(order(zi))
                Ti <- Ti[indi]
                zi <- zi[indi]
                if (ri < MAX.r[i]) {
                    if (infMAX) {
                        Ti <- rep(Ti, length.out =  MAX.r[i])
                        Ti[(ri + 1L):MAX.r[i]] <- NA
                        zi <- c(zi, rep(-Inf, MAX.r[i] - ri))
                        ri <- MAX.r[i]
                    } 
                } else {
                    Ti <- Ti[1L:MAX.r[i]]
                    zi <- zi[1L:MAX.r[i]]
                    ri <- MAX.r[i]
                }
            } else {
                ri <- 0
                if (infMAX & (0 < MAX.r[i])) {
                    ## note that here Ti has no tz attribute
                    Ti <- rep(NA, length.out = MAX.r[i])
                    zi <- rep(-Inf, MAX.r[i])
                    ri <- MAX.r[i]
                }
            }
            ## if there is at least an event: either real or fake (with -Inf level)
            if (ri > 0) { 
                if (first) {
                    MAXdata <- data.frame(block = blockNames[rep(i, ri)],
                                          date = Ti,
                                          var = zi,
                                          comment = rep("", ri))
                    first <- FALSE
                } else {
                    MAXdata <- rbind(MAXdata,
                                     data.frame(block = blockNames[rep(i, ri)],
                                                date = Ti,
                                                var = zi,
                                                comment = rep("", ri)))                               
                }
            }
            MAX.r[i] <- ri
        } else {
            MAX.r[i] <- NA
        }
        if (!is.null(OTmissing)) {
            ind2 <- (Breaks2 >= Breaks[i]) & (Breaks2 < Breaks[i + 1L])
            x <- as.numeric(interp2$x[ind2] - Breaks[i], units = "days") 
            y <- interp2$y[ind2]
            ## drop "feb, 29" in 'probMiss'
            if (length(y) == 366) y <- y[-60] 
            probMiss <- probMiss + y[1:365]
        }
    }

    names(MAXdata) <- c("block", "date", vars, "comment")

    probMiss <- probMiss / nBlock
    MAXinfo <- data.frame(MAXinfo, r = MAX.r)
    
    ##===========================================================================
    ## When 'dataFrames' is TRUE, the full list of blocks (years) is given in
    ## 'MAXinfo', and the levlels of the 'block' factor must correspond match
    ## the rows of 'MAXinfo'.
    ## When 'dataFrames' is FALSE, this is different.
    ##===========================================================================

    MAXdata[["block"]] <- factor(MAXdata[["block"]], levels = blockNames)
    
    if (dataFrames) {
        res <- list(MAXinfo = MAXinfo,
                    MAXdata = MAXdata,
                    monthGapStat = ms,
                    monthGapTS = mTS)
    } else {
        MAXdata <- as.list(MAXdata)
        ub <- as.logical(keepBlock)
        if (!infMAX) {
            ub <- ub & as.logical(table(MAXdata$block)) 
        }

        effDuration <- MAXinfo[ub, "duration"]
        r <- MAXinfo[ub, "r"]
        names(effDuration) <- names(r) <- rownames(MAXinfo[ub, ])
        ## drop unused levels
        MAXdata[["block"]] <- factor(MAXdata[["block"]])
        
        res <- list(effDuration = effDuration,
                    ## block = factor(MAXdata[["block"]]),
                    r  = r,
                    data = tapply(MAXdata[[vars]], MAXdata$block, list),
                    monthGapStat = ms,
                    monthGapTS = mTS)
        if (plot) MAXdata <- as.data.frame(MAXdata)
    }
    
    res$probMissing <- probMiss
    
    if (plot) {
        
        if (plotType == "max") {
            cols <- c("white", "lightcyan")
            y <- MAXdata[ , vars]
            ry <- range(y[is.finite(y) & !is.na(y)])
            plot(range(Breaks), ry, 
                 xlab = "Block", ylab = vars[1], type = "n", ...)
            py <- par()$usr[3:4]
            py <- py + diff(py) * c(1, -1) / 200
            for (i in 1L:nBlock) {
                rect(xleft = Breaks[i], xright = Breaks[i + 1],
                     ybottom = py[1], ytop = py[2], col = cols[1 + i %% 2], border = NA)
            }
            points(x = MAXdata[["date"]], y = MAXdata[[vars]],
                   type = "p", pch = 16, cex = 0.8)
            lines(x = MAXdata[["date"]], y = MAXdata[[vars]], type = "h")
        } else {
            if (!monthGapStat) {
                stop("the value \"gaps\" for 'plotType' works only when ",
                     "'montStat' is TRUE")
            }
                        
            days <- seq(from = as.POSIXct("2001-01-01", tz = TZ),
                        by = "day", length.out = 365)

            monthMid <- seq(from = as.POSIXct("2001-01-15", tz = TZ),
                            by = "month", length.out = 12)
            yEnd <- rep(2001, nrow(ms))
            yEnd[as.integer(ms$month) == 12] <- 2002
            dStart <- as.POSIXct(format(ms$start, format = "2001-%m-%d"),
                                 tz = TZ)
            dEnd <- as.POSIXct(paste(yEnd,
                                     format(ms$end, format = "-%m-%d"), sep = ""),
                               tz = TZ)
            
            durMiss <- with(ms, tapply(durMissing, month, sum))
            dur <- with(ms, tapply(dur, month, sum))
            monthRate <- durMiss / dur
      
            plot(days, res$probMissing, type = "n", ylim = c(0, 1),
                 xlab = "day in year", ylab = "prob. gap", xaxt = "n",
                 ...)
            abline(v = monthStart)
            mtext(side = 1, line = 0.6, at = monthMid, text = monthLab)
            set.seed(jitterSeed)
            msMissing <- ms$missing + rnorm(ms$missing, sd = 0.01)
            barCol <-  col2rgb("gray") / 256
            barCol <- rgb(red = barCol[1], green = barCol[1], blue = barCol[1],
                          alpha = 0.4)
            segments(x0 = dStart,
                     x1 = dEnd,
                     y0 = msMissing,
                     y1 = msMissing,
                     col = barCol)
            lines(x = c(monthStart[1], rep(monthStart[2:12], each = 2),
                      monthStart[13]),
                  y = rep(monthRate, each = 2),
                  col = "orange", lwd = 2) 
            lines(days, res$probMissing, type = "l",
                  lwd = 2, col = "orangered")
        }
    }
    
    res  
    
}




