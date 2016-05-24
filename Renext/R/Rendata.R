##**********************************************************************
## Authors: Yves Deville  <deville.yves@alpestat.com>
##          IRSN          <renext@irsn.fr>
##
## Functions for data managment. Some are exported, other not.
##
##**********************************************************************


##======================================================================
## NOT exported at the time
##
## read and check MAXdata from given inputs
## Note that x is unused unless it has class 'Rendata' to mimic the
## 'fitRenouv' entries
##======================================================================

makeMAXdata <- function(x,
                        data = NULL,
                        effDuration = NULL,
                        trace = 0) {
  
    flag <- FALSE
    blockNames <- NULL
    dataNames <- NULL
    
    if (missing(data)) {
        
        if ( is(x, "Rendata") && !is.null(x$MAXdata) ) {
            
            vn <- x$info$varName
            if (trace) cat("using 'MAXdata' within a 'Rendata' object\n")
            effDuration <- x$MAXinfo$duration
            block <- factor(x$MAXdata$block, levels = 1L:nrow(x$MAXinfo))
            data <- tapply(x$MAXdata[ , vn], x$MAXdata$block, list)
            ## names(MAXdata) <- paste("block", 1:length(MAXdata), sep = "")
            
            if (length(effDuration) !=  length(data))
                stop("in object 'x', elements 'MAXinfo' and 'MAXdata' do not have ",
                     "the same number of blocks")
            flag <- TRUE

            ## try to extract block names
            if ( !is.null(x$MAXinfo$comment) &&
                all(nchar(as.character(x$MAXinfo$comment)) > 0L) ){
                blockNames <- as.character(x$MAXinfo$comment)
            } else if (!any(is.na(x$MAXinfo$start)) && !any(is.na(x$MAXinfo$end))) {
                blockNames <- paste(formatPeriod(start = x$MAXinfo$start,
                                                 end = x$MAXinfo$end), "(MAX)")
            } else blockNames <- paste("MAX block", 1L:nrow(x$MAXinfo))

            ## try to extract data names
            if ( !is.null(x$MAXdata$comment) &&
                any(nchar(as.character(x$MAXdata$comment)) > 0L) ){
                dataNames <- as.character(x$MAXdata$comment)
            } else if (!any(is.na(x$MAXdata$date))) {
                dataNames <- format(x$MAXdata$date, "%Y-%m-%d")
            } else dataNames <- rep("", nrow(x$MAXdata))
            
            
        } else {
            flag <- FALSE
            block <- NULL
            effDuration <- NULL
            threshold <- NULL
            r <- NULL
        }
        
    } else {
        
        if (is.null(data)) {
            flag <- FALSE
            block <- NULL
            effDuration <- NULL
            threshold <- NULL;
            r <- NULL
        } else if (is.list(data)) {
            
            if (length(data) != length(effDuration))
                stop("when 'data' is a list, 'effDuration' must be a numeric ",
                     "vector of length length(data)")

            blockNames <- names(data)
            if (is.null(blockNames)) blockNames <- paste("MAX block", 1L:length(data))
            
            r <- unlist(lapply(data, length))
            block <- factor(rep(1:length(data), times = r), levels = 1:length(data))
            dataNames <- lapply(data, fillNames)
            
            flag <- TRUE
            
        } else {
            if (length(effDuration) != 1L)
                stop("when 'data' is a vector, 'effDuration' must be a numeric ",
                     "with length 1")
            block <- factor(rep(1L, times = length(data)), levels = 1)
            r <- length(data)
            dataNames <- fillNames(data)
            blockNames <- "MAX block1" 
            data <- list("MAX block 1" = data)
            flag <- TRUE
        }
        
    }
    
    if (!flag) {
        block <- NULL
        effDuration <- NULL
        r <- NULL
    } else {
        r <- unlist(lapply(data, length))
        names(effDuration) <- names(r) <- paste("block", 1L:length(effDuration))
    }
    
    list(flag = flag,
         block = block,
         blockNames = blockNames,
         effDuration = effDuration,
         r = r,
         data = data,
         dataNames = dataNames)    
    
}

##======================================================================
## read and check OTSdata from given inputs
## Note that x is unused unless it has class 'Rendata'
##======================================================================

makeOTSdata <- function(x,
                        data = NULL,
                        effDuration = NULL,
                        threshold = NULL,
                        trace = 0) {
    
    flag <- FALSE
    blockNames <- NULL
    dataNames <- NULL
    
    if (missing(data)) {
        
        if ( is(x, "Rendata") && !is.null(x$OTSdata) ) {
            
            vn <- x$info$varName
            if (trace) cat("using 'OTSdata' within a 'Rendata' object\n")
            effDuration <- x$OTSinfo$duration
            threshold <- x$OTSinfo$threshold
            r <- x$OTSinfo$r
            
            block <- factor(x$OTSdata$block, levels = 1:nrow(x$OTSinfo))
            data <- tapply(x$OTSdata[ , vn], block, as.numeric)
            ## to get empty items as numeric(0) rather than NULL...
            data <- lapply(data, function(x) as.numeric(x))
            ind <- (table(block) == 0)
            if (any(ind)) {
                for (i in (1L:length(data))[ind]) data[[i]] <- numeric(0)
            }
            
            flag <- TRUE
            
            ## try to extract block names
            if ( !is.null(x$OTSinfo$comment) &&
                all(nchar(as.character(x$OTSinfo$comment)) > 0L) ){
                blockNames <- as.character(x$OTSinfo$comment)
            } else if (!any(is.na(x$OTSinfo$start)) && !any(is.na(x$OTSinfo$end))) {
                blockNames <- paste(formatPeriod(start = x$OTSinfo$start,
                                                 end = x$OTSinfo$end), "(OTS)")
            } else blockNames <- paste("OTS block", 1L:nrow(x$OTSinfo))

            ## try to extract data names
            if ( !is.null(x$OTSdata$comment) &&
                any(nchar(as.character(x$OTSdata$comment)) > 0L) ){
                dataNames <- as.character(x$OTSdata$comment)
            } else if (!any(is.na(x$OTSdata$date))) {
                dataNames <- format(x$OTSdata$date, "%Y-%m-%d")
            } else dataNames <- rep("", nrow(x$OTSdata))
            
        } else {
            if (trace) cat("'x' not of class 'Rendata': ignored\n")
            flag <- FALSE
        }
    } else {
        
        if (is.null(data)) {
            flag <- FALSE
            block <- NULL;
            effDuration <- NULL
            threshold <- NULL;
            r <- NULL
        } else if (is.list(data)) {
            if ( (length(data) != length(effDuration)) ||
                (length(data) != length(threshold)) )
                stop("when 'data' is a list, 'effDuration' and 'threshold' must be a ",
                     "numeric vector of length length(data)")
            blockNames <- names(data)
            if (is.null(blockNames)) blockNames <- paste("OTS block", 1L:length(data))
            r <- unlist(lapply(data, length))
            block <- factor(rep(1L:length(data), times = r), levels = 1L:length(data))
            dataNames <- lapply(data, fillNames)
            flag <- TRUE
        } else {
            if ((length(effDuration) != 1) || (length(threshold) != 1L))
                stop("when 'data' is a vector, 'effDuration' and 'threshold' must be ",
                     "numeric with length 1")
            r <- length(data)
            block <- factor(1L, levels = 1)
            dataNames <- fillNames(data)
            blockNames <- "OTS block1" 
            data <- list("OTS block1" = data)
            flag <- TRUE
        }
        
    }
    
    if (!flag) {
        block <- NULL
        effDuration <- NULL
        threshold <- NULL
        r <- NULL
    } else {
        names(effDuration) <- names(r) <-  names(threshold) <-
            paste("block", 1:length(effDuration))
    }
    
    list(flag = flag,
         block = block,
         blockNames = blockNames,
         effDuration = effDuration,
         threshold = threshold,
         r = r,
         data = data,
         dataNames = dataNames)    
    
}

##============================================================
## plot composite datasets as those read by readXML function
## from an XML/csv source.
##
##
##============================================================

plot.Rendata <- function(x,
                         textOver = quantile(x$OTdata[, x$info$varName], probs = 0.99),
                         showHist = TRUE,
                         ...) {
    
    ## this is only to avoid a NOTE in build step
    ## " no visible binding for global variable 'block' "
    block <- 1 
    
    periodsLeg <- c("OTdata", "OTSdata", "MAXdata")
    periodsBg <-  c("lightcyan", "lightyellow", "DarkOliveGreen1")
    names(periodsBg) <- periodsLeg
    periodsFg <- c("cyan", "gold", "DarkOliveGreen2")
    names(periodsFg) <- periodsLeg
    periodsFlag <- c(TRUE, FALSE, FALSE)
    names(periodsFlag) <- periodsLeg
    
    y <- x$OTdata[ , x$info$varName]
    start <- as.POSIXct(x$OTinfo$start)
    end <- as.POSIXct(x$OTinfo$end)
    ymin <-  x$OTinfo$threshold ## min(y)
    ymax <- max(y)
    
    if (showHist) {
        if (!is.null(x$OTSinfo)) {
            OTSstart <- as.POSIXct(x$OTSinfo$start)
            OTSend <- as.POSIXct(x$OTSinfo$start)
            start <- min(start, OTSstart)
            end <- max(end, OTSend)
            if (!is.null(x$OTSdata)) {
                ymin <- min(ymin, min(x$OTSdata[ , x$info$varName]))
                ymax <- max(ymax, max(x$OTSdata[ , x$info$varName]))
            }
        }
        
        if (!is.null(x$MAXinfo)) {
            MAXstart <- as.POSIXct(x$MAXinfo$start)
            MAXend <- as.POSIXct(x$MAXinfo$start)
            start <- min(start, MAXstart)
            end <- max(end, MAXend)
            ymin <- min(ymin, min(x$MAXdata[ , x$info$varName]))
            ymax <- max(ymax, max(x$MAXdata[ , x$info$varName]))
        }
    }
    
    yLab <- x$info$varName
    if (!is.null(x$info$varUnit))
        yLab <- paste(yLab, " (", x$info$varUnit, ")", sep = "")
    
    plot(x = x$OTdata[ , "date"],
         y = x$OTdata[ , x$info$varName],
         type ="n",
         xlim = c(start, end),
         ylim = c(ymin, ymax),
         xlab = " ",
         ylab = yLab,
         main = x$info$longLab,
         ...)
    
    rg <- par()$usr[3:4]
    drg <- rg[2]-rg[1]
    rg <- rg + c(drg, -drg)/100
    ## cat("rg = \n"); print(rg)
    
    rect(xleft = x$OTinfo$start,
         xright = x$OTinfo$end,
         ybottom = rg[1],
         ytop = rg[2],
         col = periodsBg["OTdata"], border = periodsFg["OTdata"])
    
    ## show the missing periods
    if (!is.null(x$OTmissing)) {
        for (i in 1:nrow(x$OTmissing)) {
            polygon(x = c(x$OTmissing$start[i], x$OTmissing$end[i],
                        x$OTmissing$end[i],  x$OTmissing$start[i]),
                    y = rep(rg, each = 2),
                    border = periodsBg["OTdata"], col = "white")
        }
    }
    
    lines(x = x$OTdata[ , "date"],
          y = x$OTdata[ , x$info$varName],
          type ="h",
          col = "SteelBlue3")
    
    ## show the threshold for OTdata
    lines(x = as.POSIXct(c(x$OTinfo$start, x$OTinfo$end)),
          y = rep(as.numeric(x$OTinfo$threshold), times = 2),
          col = "orange", lwd = 2)
    
    if (!is.na(textOver)){
        ind <- (x$OTdata[, x$info$varName] > textOver)
        if (any(ind)) {
            points(x = x$OTdata[ind, "date"],
                   y = x$OTdata[ind, x$info$varName],
                   pch = 21, cex = 0.7,
                   col = "SteelBlue4")
            text(x = x$OTdata[ind, "date"],
                 y = x$OTdata[ind, x$info$varName],
                 labels = format(x$OTdata[ind, "date"], "%Y-%m-%d"),
                 col = "SteelBlue4",
                 pos = 4, cex = 0.7)
        }
    }
    
    ##=================================================================
    ## plot OTS data  if present
    ##=================================================================
    
    if (showHist && !is.null(x$OTSinfo)) {
        periodsFlag["OTSdata"] <- TRUE
        
        ## for each 'OTS' block ..; 
        for (i in 1:nrow(x$OTSinfo)) {
            
            rect(xleft = x$OTSinfo$start[i],
                 xright = x$OTSinfo$end[i],
                 ybottom = rg[1],
                 ytop = rg[2],
                 col = periodsBg["OTSdata"], border = periodsFg["OTSdata"])
            
            lines(x = as.POSIXct(c(x$OTSinfo$start[i], x$OTSinfo$end[i])),
                  y = rep(x$OTSinfo$threshold[i], times = 2),
                  col = "orange", lwd = 2)
            
        }
        
        if (!is.null(x$OTSdata)) {
            
            for (i in 1:nrow(x$OTSinfo)) {
                
                datai <- subset(x$OTSdata, block == i) 
                ind <- !is.na(datai[ , "date"])
                
                if (any(ind)) {
                    lines(x = datai[ind, "date"],
                          y = datai[ind, x$info$varName],
                          type ="h",
                          col = "red3")
                }
                if (any(!ind)) {
                    segments(x0 = rep(x$OTSinfo$start[i], sum(!ind)),
                             x1 = rep(x$OTSinfo$end[i], sum(!ind)),
                             y0 = datai[!ind , x$info$varName],
                             y1 = datai[!ind , x$info$varName],
                             lty = "dashed",
                             col = "red3")
                }
                
            }
            
            if (!is.na(textOver)){
                ## Added checks for NA on 2013-06-14
                ## ind <- (x$OTSdata[, x$info$varName] > textOver)
                ind <- (x$OTSdata[ , x$info$varName] > textOver) &
                    !is.na(x$OTSdata[ , "date"])
                if (any(ind)) {
                    points(x = x$OTSdata[ind, "date"],
                           y = x$OTSdata[ind, x$info$varName],
                           pch = 21, cex = 0.7,
                           col = "red3")
                    text(x = x$OTSdata[ind, "date"],
                         y = x$OTSdata[ind, x$info$varName],
                         labels = format(x$OTSdata[ind, "date"], "%Y-%m-%d"),
                         col = "red3",
                         pos = 4, cex = 0.7)
                }
            }
        }
    }
    
    ##=================================================================
    ## plot MAX data if present
    ## Note that 'MAX' historical data blocks alway contain observation,
    ## which is not true for 'OTS' data.
    ##=================================================================
    
    if (showHist && !is.null(x$MAXinfo)) {
        
        periodsFlag["MAXdata"] <- TRUE
        
        if (!is.null(x$MAXdata)) {
            
            ## for each 'MAX' block ..; 
            for (i in 1:nrow(x$MAXinfo)) {
                rect(xleft = x$MAXinfo$start[i],
                     xright = x$MAXinfo$end[i],
                     ybottom = rg[1],
                     ytop = rg[2],
                     col = periodsBg["MAXdata"], border = periodsFg["MAXdata"])  
                
                datai <- subset(x$MAXdata, block == i) 
                ind <- !is.na(datai[ , "date"])
                
                if (any(ind)) {
                    lines(x = datai[ind, "date"],
                          y = datai[ind, x$info$varName],
                          type ="h",
                          col = "SpringGreen4")
                }
                if (any(!ind)) {
                    segments(x0 = rep(x$MAXinfo$start[i], sum(!ind)),
                             x1 = rep(x$MAXinfo$end[i], sum(!ind)),
                             y0 = datai[!ind , x$info$varName],
                             y1 = datai[!ind , x$info$varName],
                             lty = "dashed",
                             col = "SpringGreen4")
                }
            }
            
            if (!is.na(textOver)){
                
                ind <- (x$MAXdata[ , x$info$varName] > textOver) &
                    !is.na(x$MAXdata[ , "date"])
                ## nothing shown for missing dates... 
                if (any(ind)) {
                    points(x = x$MAXdata[ind, "date"],
                           y = x$MAXdata[ind, x$info$varName],
                           pch = 21, cex = 0.7,
                           col = "SpringGreen4")
                    text(x = x$MAXdata[ind, "date"],
                         y = x$MAXdata[ind, x$info$varName],
                         labels = format(x$MAXdata[ind, "date"], "%Y-%m-%d"),
                         col = "SpringGreen4",
                         pos = 4, cex = 0.7)
                }
            }
        }
    }
    
    legend("topleft",
           fill = periodsBg[periodsFlag],
           col = periodsFg[periodsFlag],
           legend = periodsLeg[periodsFlag])
    
}

##========================================================================
## summary method for class 'Rendata'
##
## We keep information from 'info', 'OTinfo' and summarize the 'data'
## parts
##
##========================================================================

summary.Rendata <- function(object, ...) {

    ans <- list(info = object$info,
                OTinfo = object$OTinfo,
                MAXinfo = object$MAXinfo,
                OTSinfo = object$OTSinfo)
    
    ans$info <- paste(sprintf("o Dataset %s", object$info$shortLab),
                      sprintf("   data '%s', variable '%s' (%s)",
                              object$info$name, object$info$varName, object$info$varUnit),
                      sep = "\n")  
    
    ans$OTinfo <- paste("o OT data (main sample)",
                        "from ", format(object$OTinfo$start, "%Y-%m-%d"),
                        " to ", format(object$OTinfo$end, "%Y-%m-%d"),
                        sprintf(" (eff. dur. %6.2f years)\n", object$OTinfo$effDuration))
    
    var <- object$OTdata[ , object$info$varName]
    ans$OTsummary <- c(n = length(var), summary(var))
    
    if (!is.null(object$OTmissing)) {
        dur <- as.numeric(as.POSIXct(object$OTmissing$end) -
                              as.POSIXct(object$OTmissing$start),
                          units = "days")
        ans$OTmissing <- sprintf("o missing 'OT' periods, total %5.2f years",
                                 sum(dur)/365.25)
        ans$OTmissingsummary <- c(n = length(dur), summary(dur/365.25))
    } else   {
        ans$OTmissing <- "o no missing OT periods"
        ans$OTmissingsummary <- NULL
    }
    
    if (!is.null(object$MAXinfo)) {
        ans$MAXinfo <-
            sprintf("o 'MAX' historical info: %d blocks, %d obs., total duration = %5.2f years",
                    nrow(object$MAXinfo), nrow(object$MAXdata), sum(object$MAXinfo$duration))
    } else   ans$MAXinfo <- "o no 'MAX' historical data"
    
    
    if (!is.null(object$OTSinfo)) {
        ans$OTSinfo <-
            sprintf("o 'OTS' historical info: %d blocks, %d obs., total duration = %5.2f years",
                    nrow(object$OTSinfo), nrow(object$OTSdata), sum(object$OTSinfo$duration))
    } else   ans$OTSinfo <- "o no 'OTS' historical data"
    
    
    class(ans) <- "summary.Rendata"
    ans
    
}

print.summary.Rendata <- function(x, ...) {
    cat(x$info, "\n")
    cat("\n")
    cat(x$OTinfo, "\n")
    print(x$OTsummary)
    cat("\n")
    cat(x$OTmissing, "\n\n")
    if (!is.null(x$OTmissingsummary)){
        print(x$OTmissingsummary)
        cat("\n")
    }
    cat(x$MAXinfo, "\n\n")
    cat(x$OTSinfo, "\n\n")
}

print.Rendata <- function(x, ...) {
    print(summary(x, ...))
}
