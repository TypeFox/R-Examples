"plot.drc" <-
function(x, ..., add = FALSE, level = NULL, type = c("average", "all", "bars", "none", "obs", "confidence"), 
broken = FALSE, bp, bcontrol = NULL, conName = NULL, axes = TRUE, gridsize = 100, 
log = "x", xtsty, xttrim = TRUE, xt = NULL, xtlab = NULL, xlab, xlim, 
yt = NULL, ytlab = NULL, ylab, ylim,
cex, cex.axis = 1, col = FALSE, lty, pch, 
legend, legendText, legendPos, cex.legend = 1,
normal = FALSE, normRef = 1, confidence.level = 0.95)
{    
    object <- x
    type <- match.arg(type)
    
    ## Determining logarithmic scales
    if ((log == "") || (log == "y"))
    {
        logX <- FALSE
    } else {
        logX <- TRUE
    }    
    
    ## Determining the tick mark style for the dose axis
    if (missing(xtsty))
    {
        if (logX)
        {
            xtsty <- "base10"
        } else {
            xtsty <- "standard"
        }
    }

    ## Constructing the plot data
#    origData<-object$"data"
#    doseDim <- ncol(origData) - 4  # subtracting 4 because the data frame also contains columns with response, 
                                   # curve no. (in new and old enumeration) and weights
#    if (doseDim > 1) {stop("No plot features for plots in more than two dimensions")}
#    dose1 <- origData[, 1:doseDim]
#    resp1 <- origData[, doseDim+1]
    
    dataList <- object[["dataList"]]
    dose <- dataList[["dose"]]
    resp <- dataList[["origResp"]]
    curveid <- dataList[["curveid"]]
    plotid <- dataList[["plotid"]]    
#    levels(curveid) <- unique(dataList[["curveid"]])  # commented out 12/2-2014
    
    ## Normalizing the response values
    if (normal)
    {
        respList <- split(resp, curveid)
        resp <- unlist(mapply(normalizeLU, respList, as.list(as.data.frame(getLU(object))), normRef = normRef))            
#        respNew <- unlist(mapply(normalizeLU, respList, as.list(as.data.frame(getLU(object)))))    
#        print(respNew)
#        print(resp)
    }
    
#    assayNo <- origData[, 3]
#    assayNoOld <- as.vector(origData[, 4])  # as.vector() to remove factor structure

#    if (!is.null(dataList[["plotid"]]))
#    {
#        assayNoOld <- dataList[["plotid"]]
#        uniAss <- unique(assayNoOld)        
#    } else {
#        assayNoOld <- as.vector(curveid)  # as.vector() to remove factor structure
#        uniAss <- unique(curveid)
#    }   
# commented out 3/3-'14 
    
    if (!is.null(plotid))           
    {  # used for event-time data
    
#        piVec <- dataList[["plotid"]]
#        levels(piVec) <- as.vector(unique(curveid))
#        assayNoOld <- as.vector(piVec)  # used for event-time data         
        assayNoOld <- as.vector(plotid)  
    } else {
        assayNoOld <- as.vector(curveid)
    }
    uniAss <- unique(assayNoOld) 
    numAss <- length(uniAss) 

#    print(uniAss)
    
#    numAss <- length(unique(assayNoOld))
#    doPlot <- is.null(level) || any(unique(assayNoOld) %in% level)

    doPlot <- is.null(level) || any(uniAss %in% level)
    if (!doPlot) {stop("Nothing to plot")}

    plotFct <- (object$"curve")[[1]]
    logDose <- (object$"curve")[[2]]
    naPlot <- ifelse(is.null(object$"curve"$"naPlot"), FALSE, TRUE)

#    if (missing(conName)) 
#    {
#        if (is.null(logDose)) {conName <- expression(0)} else {conName <- expression(-infinity)}
#    }

    ## Assigning axis names
#    varNames <- colnames(origData)[1:(doseDim+1)]  

    dlNames <- dataList[["names"]]
    doseName <- dlNames[["dName"]]
    respName <- dlNames[["orName"]]
    # axis names are the names of the dose variable and response variable in the original data set
#    if (missing(xlab)) {if (varNames[1] == "") {xlab <- "Dose"} else {xlab <- varNames[1]}}
#    if (missing(ylab)) {if (varNames[2] == "") {ylab <- "Response"} else {ylab <- varNames[2]}}
    if (missing(xlab)) {if (doseName == "") {xlab <- "Dose"} else {xlab <- doseName}}
    if (missing(ylab)) {if (respName == "") {ylab <- "Response"} else {ylab <- respName}}     

    ## Determining range of dose values
    if (missing(xlim)) 
    {
        xLimits <- c(min(dose), max(dose))
    } else {
        xLimits <- xlim  # if (abs(xLimits[1])<zeroEps) {xLimits[1] <- xLimits[1] + zeroEps}
    }

    ## Handling small dose values
    if (missing(bp)) 
    {
#        conLevel <- ifelse(is.null(logDose), 1e-2, log(1e-2))
        
        ## Constructing appropriate break on dose axis
        if (!is.null(logDose))  # natural logarithm
        {
            conLevel <- round(min(dose[is.finite(dose)])) - 1
        } else {
            log10cl <- round(log10(min(dose[dose > 0]))) - 1
            conLevel <- 10^(log10cl)
        }
    } else {
        conLevel <- bp
    }   
    if ((xLimits[1] < conLevel) && (logX || (!is.null(logDose)))) 
    {
        xLimits[1] <- conLevel
        smallDoses <- (dose < conLevel)
        dose[smallDoses] <- conLevel
        if (is.null(conName)) 
        { 
            if (is.null(logDose)) {conName <- expression(0)} else {conName <- expression(-infinity)}
        }
#        conNameYes <- TRUE 
    } else {
#        conNameYes <- FALSE
        conName <- NULL
    }
    if (xLimits[1] >= xLimits[2]) {stop("Argument 'conLevel' is set too high")}

    ## Constructing dose values for plotting
#    if (doseDim == 1) 
#    {
    if ((is.null(logDose)) && (logX))
    {
       dosePts <- exp(seq(log(xLimits[1]), log(xLimits[2]), length = gridsize))
       ## Avoiding that slight imprecision produces dose values outside the dose range
       ## (the model-robust predict method is sensitive to such deviations!)
       dosePts[1] <- xLimits[1]
       dosePts[gridsize] <- xLimits[2]           
    } else {
       dosePts <- seq(xLimits[1], xLimits[2], length = gridsize)
    }
#    } else {}  # No handling of multi-dimensional dose values


    ## Finding minimum and maximum on response scale
    if (is.null(logDose)) 
    {
        plotMat <- plotFct(dosePts)
    } else {
        plotMat <- plotFct(logDose^(dosePts))
    }
    ## Normalizing the fitted values
    if (normal)
    {
        respList <- split(resp, curveid)
        plotMat <- mapply(normalizeLU, as.list(as.data.frame(plotMat)), as.list(as.data.frame(getLU(object))), normRef = normRef)            
#        pmNew <- mapply(normalizeLU, as.list(as.data.frame(plotMat)), as.list(as.data.frame(getLU(object))))    
#        print(pmNew)
#        print(plotMat)
    }    
#    numCol <- ncol(plotMat)

    maxR <- max(resp)
    options(warn = -1)  # suppressing warning in case maximum of NULL is taken
    maxPM <- apply(plotMat, 2, max, na.rm = TRUE)
    if (max(maxPM) > maxR) {maxPM <- maxPM[which.max(maxPM)]} else {maxPM <- maxR}
    options(warn=0)  

    if (missing(ylim)) 
    {
        if (missing(xlim))
        {
            yLimits <- c(min(resp), maxPM)
        } else {
            yLimits <- getRange(dose, resp, xLimits)
        }            
    } else {
        yLimits <- ylim
    }


    ## Cutting away y values (determined by the fitted model) outside the limits
###    naPlot <- FALSE  # remove naPlot further down
#    for (i in 1:numCol)
#    {
#        logVec <- !(plotMat[, i] >= yLimits[1]  & plotMat[, i] <= yLimits[2])
#        if ( any(!is.na(logVec)) && any(logVec) )
#        {
#            plotMat[logVec, i] <- NA
#            naPlot <- TRUE
#        }
#    }

    ## Setting a few graphical parameters
    par(las = 1)
    if (!is.null(logDose)) 
    {
        if (log == "x") {log <- ""}
        if ( (log == "xy") || (log == "yx") ) {log <- "y"}
    }  
    
    ## Cutting away original x values outside the limits
    eps1 <- 1e-8
    logVec <- !( (dose < xLimits[1] - eps1) | (dose > xLimits[2] + eps1) )
    dose <- dose[logVec]
    resp <- resp[logVec]
#    assayNo <- assayNo[logVec]
    assayNoOld <- assayNoOld[logVec]    
     
    ## Calculating predicted values for error bars
    if (identical(type, "bars"))
    {
#        predictMat <- predict(object, interval = "confidence")[, 3:4]
#        predictMat <- predict(object, interval = "confidence")[, c("Lower", "Upper")]
        predictMat <- predict(object, interval = "confidence",
                              level = confidence.level)[, c("Lower", "Upper")]        
#        print(predictMat)
    
        barFct <- function(plotPoints)
        {
            pp3 <- plotPoints[, 3]
            pp4 <- plotPoints[, 4]
            plotCI(plotPoints[, 1], pp3 + 0.5 * (pp4 - pp3), 
            li = pp3, ui = pp4, add = TRUE, pch = NA)
        }

        ciFct <- function(level, ...){invisible(NULL)}
        
        pointFct <- function(plotPoints, cexVal, colVal, pchVal, ...){invisible(NULL)} 

    } else if (identical(type, "confidence"))
    {
      
        barFct <- function(plotPoints){invisible(NULL)}
      
        ciFct <- function(level, ...)
        {
            newdata <- data.frame(DOSE=dosePts, CURVE=rep(level, length(dosePts)))
            predictMat <- predict(object, 
                                  newdata=newdata,
                                  interval = "confidence",
                                  level=confidence.level)
        
            x <- c(dosePts, rev(dosePts))
            y <- c(predictMat[,"Lower"], rev(predictMat[,"Upper"]))
            polygon(x,y, border=NA, ...)
        }
      
        pointFct <- function(plotPoints, cexVal, colVal, pchVal, ...){invisible(NULL)} 
        
    } else {
      
        barFct <- function(plotPoints){invisible(NULL)}
  
        ciFct <- function(level, ...){invisible(NULL)}
  
        pointFct <- function(plotPoints, cexVal, colVal, pchVal, ...)
        {
            points(plotPoints, cex = cexVal, col = colVal, pch = pchVal, ...)                
        }
    }                


    ## Setting the plot type
    if ( (identical(type, "none")) || (identical(type, "bars")) )
    {
        plotType <- "n"
    } else {
        plotType <- "p"
    }

    ## Determining levels to be plotted
#    uniAss <- unique(assayNoOld)
    if (is.null(level)) 
    {
        level <- uniAss
    } else {
        level <- intersect(level, uniAss)
    }
    lenlev <- length(level)
    
    ## Determining presence of legend
    if (missing(legend))
    {
        if (lenlev == 1) {legend <- FALSE} else {legend <- TRUE}
    }

    ## Setting graphical parameters
    colourVec <- rep(1, lenlev)
    if (is.logical(col) && col) 
    {
        colourVec <- 1:lenlev
    }  
    if (!is.logical(col) && (length(col) == lenlev) )
    {
        colourVec <- col
    }
    if (!is.logical(col) && (!(length(col) == lenlev)) ) 
    {
        colourVec <- rep(col, lenlev)
    }   
    cexVec <- parFct(cex, lenlev, 1)
    ltyVec <- parFct(lty, lenlev)
    pchVec <- parFct(pch, lenlev)           
      
    ## Plotting data
    levelInd <- 1:lenlev
    for (i in levelInd)
    {
        indVec <- level[i] == assayNoOld
#        print(indVec)
#        print(level[i])
#        print(assayNoOld)
        plotPoints <- 
        switch(
            type, 
            
            "average" = cbind(as.numeric(names(tapVec <- tapply(resp[indVec], 
            dose[indVec], mean))), tapVec),
            
            "bars"    = cbind(
            as.numeric(names(tapVec <- tapply(resp[indVec], dose[indVec], mean))), 
            tapVec, 
            tapply(predictMat[indVec, 1], dose[indVec], head, 1),
            tapply(predictMat[indVec, 2], dose[indVec], head, 1)),
            
            "none"    = cbind(dose[indVec], resp[indVec]),
            "all"     = cbind(dose[indVec], resp[indVec]),
            "obs"     = cbind(dose[indVec], resp[indVec])
        )
#        print(plotPoints)
        ## New approach
#        plotPointsRaw <- ppList[uniAss2[i]]
#        plotPoints <- with(plotPointsRaw,
#        switch(
#            type, 
#            "average" = cbind(as.numeric(names(tapVec <- tapply(resp, 
#            dose, mean))), tapVec),
#            "bars"    = cbind(as.numeric(names(tapVec <- tapply(resp, 
#            dose, mean))), tapVec, tapply(predictMat[indVec, 1], dose[indVec], head, 1),
#            tapply(predictMat[indVec, 2], dose[indVec], head, 1)),
#            "none"    = cbind(dose[indVec], resp[indVec]),
#            "all"     = cbind(dose[indVec], resp[indVec]),
#            "obs"     = cbind(dose[indVec], resp[indVec])
#        ))     
        
        if ( (!add) && (i == 1) )
        {
            ## Plotting data for the first curve id
            plot(plotPoints, type = plotType, xlab = xlab, ylab = ylab, log = log, xlim = xLimits, ylim = yLimits, 
            axes = FALSE, frame.plot = TRUE, col = colourVec[i], pch = pchVec[i], cex = cexVec[i], ...) 
            
            ## Adding error bars
            barFct(plotPoints)      
            
            ## Add confidence regions
            ciFct(level=i, col=alpha(colourVec[i],0.25))            
            
            ## Adding axes    
            addAxes(axes, cex.axis, conName, xt, xtlab, xtsty, xttrim, logX, yt, ytlab, conLevel, logDose)

            ## Adding axis break
            ivMid <- brokenAxis(bcontrol, broken, conLevel, dosePts, gridsize, log, logX, logDose)       
                           
        ## Plotting in the case "add = TRUE" and for all remaining curve ids
        } else {  
            ## Adding axis break (in fact only restricting the dose range to be plotted)
            ivMid <- brokenAxis(bcontrol, broken, conLevel, dosePts, gridsize, log, logX, logDose, plotit = FALSE)       
            
            if (!identical(type, "none"))  # equivalent of type = "n" in the above "plot" 
            {
                pointFct(plotPoints, cexVec[i], colourVec[i], pchVec[i], ...)
            
                ## Adding error bars
                barFct(plotPoints)
                
                ## Add confidence regions
                ciFct(level=i, col=alpha(colourVec[i],0.25))
            }
        }
    }
    
    ## Plotting fitted curves
    noPlot <- rep(FALSE, lenlev)
    if (!identical(type, "obs"))
    {
        for (i in levelInd)
        {
            indVal <- uniAss %in% level[i]
            if ( (!naPlot) && (any(is.na(plotMat[, indVal]))) ) 
            {
                noPlot[i] <- TRUE
                next
            }
            lines(dosePts[ivMid], plotMat[ivMid, indVal], lty = ltyVec[i], col = colourVec[i], ...)
        }
    }
    

#    ## Plotting pointwise prediction intervals
#    if (identical(type, "predict"))
#    {
#        for (i in levelInd)
#        {
#            indVal <- uniAss %in% level[i]
#            if ( (!naPlot) && (any(is.na(plotMat[, indVal]))) ) 
#            {
#                noPlot[i] <- TRUE
#                next
#            }
#            lines(dosePts[ivMid], plotMat[ivMid, indVal], lty = ltyVec[i], col = colourVec[i], ...)
#        }
#    }

    
    ## Adding legend
    makeLegend(colourVec, legend, cex.legend, legendPos, legendText, lenlev, level, ltyVec, 
    noPlot, pchVec, type, xLimits, yLimits)
    
    ## Resetting graphical parameter
    par(las = 0) 

    retData <- data.frame(dosePts, as.data.frame(plotMat))
#    colnames(retData) <- c(colnames(origData)[1:doseDim], as.character(unique(assayNoOld)))
    colnames(retData) <- c(doseName, dlNames[["cNames"]])
    
    invisible(retData)
}


"getRange" <- function(x, y, xlim)
{
    logVec <- ((x >= xlim[1]) & (x <= xlim[2]))
    return(range(y[logVec]))
}

# if (FALSE)
# {
# "addAxes" <- function(axes, cex.axis, conLevel, conName, conNameYes, xt, xtlab, xsty, yt, ytlab)
# {
#     if (!axes) {return(invisible(NULL))}  # doing nothing
#     
#     ## Concerning the y axis
#     yaxisTicks <- axTicks(2)
#     yLabels <- TRUE
#     if (!is.null(yt)) {yaxisTicks <- yt; yLabels <- yt}
#     if (!is.null(ytlab)) {yLabels <- ytlab}                
#                 
#     ## Concerning the x axis                
#     xaxisTicks <- axTicks(1)
#     
#     
#     
#     
#     if (conNameYes)
#     {
#         xaxisTicks[1] <- conName
#     }
#     xaxisTicksOrig <- xaxisTicks
#     xLabels <- xaxisTicks
#     
#     if ((conNameYes) && (min(xt) < conLevel)) 
#     {
#         xLimits[1] <- conLevel
#         smallDoses <- dose<conLevel
#         dose[smallDoses] <- conLevel
#         conNameYes <- TRUE 
#     }
#     
# #    print(xaxisTicks)
# 
# #    ## Avoiding too many tick marks on the x axis
# #    lenXT <- length(xaxisTicks)
# #    if (lenXT > 6) 
# #    {
# #        halfLXT <- floor(lenXT/2) - 1
# #        chosenInd <- 1 + 2*(0:halfLXT)  # ensuring that control always is present
# #        xaxisTicks <- xaxisTicks[chosenInd]
# #        xLabels <- xLabels[chosenInd]
# #    }    
# #    
#     if (identical(xsty, "base10"))
#     {    
#         ceilingxTicks <- ceiling(log10(xaxisTicks[-1]))
#         xaxisTicks <- c(xaxisTicks[1], 10^(unique(ceilingxTicks)))
#         xLabels <- c(xLabels[1], unlist(tapply(xLabels[-1], ceilingxTicks, head, 1)))
# 
#         ## Reverting to original tick marks in case too few were created
#         if (length(xaxisTicks) < 3)
#         { 
#             xaxisTicks <- xaxisTicksOrig
#         }
#     }
#                 
#     if (!is.null(xt)) 
#     {
#         if (as.character(xt[1]) == as.character(eval(conName))) 
#         {
#              xaxisTicks <- c(xaxisTicks[1], xt[-1])
#              xLabels <- c(conName, xt[-1])                        
#         } else {
#              xaxisTicks <- xt
#              xLabels <- xt
#              conNameYes <- FALSE
#         }
#     }
#     if (!is.null(xtlab)) {xLabels <- xtlab}
# #    print(xaxisTicks)
# 
#     ## Updating x axis labels
#     xLabels <- as.expression(xaxisTicks)
#     if (conNameYes) {xLabels[1] <- conName}                                                
# 
#     axis(1, at = xaxisTicks, labels = xLabels, cex.axis = cex.axis)
#     axis(2, at = yaxisTicks, labels = yLabels, cex.axis = cex.axis)        
# }
# }


"addAxes" <- function(axes, cex.axis, conName, xt, xtlab, xtsty, xttrim, logX, yt, ytlab, conLevel, logDose)
{
    if (!axes) {return(invisible(NULL))}  # doing nothing
    
    ## Setting up the y axis tick mark locations and labels 
    yaxisTicks <- axTicks(2)
    yLabels <- TRUE
    if (!is.null(yt)) {yaxisTicks <- yt; yLabels <- yt}
    if (!is.null(ytlab)) {yLabels <- ytlab}                
                
    ## Setting up the x axis tick mark locations and labels               
    if (!is.null(xt)) 
    {
        xaxisTicks <- xt
        if (identical(as.numeric(xaxisTicks)[1], 0))
        { 
            xaxisTicks[1] <- conLevel
        }
    } else {
        xaxisTicks <- axTicks(1)

        ## Styling the x axis tick marks
        if (identical(xtsty, "base10"))
        {    
            if (!is.null(logDose)) 
            {
                ceilingxTicks <- ceiling(xaxisTicks[-1])
                xaxisTicksOrig <- xaxisTicks
                xaxisTicks <- c(xaxisTicks[1], unique(ceilingxTicks))
            } else {
                ceilingxTicks <- ceiling(log10(xaxisTicks[-1]))
                xaxisTicksOrig <- xaxisTicks
                xaxisTicks <- c(xaxisTicks[1], 10^(unique(ceilingxTicks)))
#        xLabels <- c(xLabels[1], unlist(tapply(xLabels[-1], ceilingxTicks, head, 1)))
            }

            ## Reverting to original tick marks in case too few were created
            if (length(xaxisTicks) < 3)
            { 
                xaxisTicks <- xaxisTicksOrig
#                xLabels <- as.character(xaxisTicks)
            }
        }
    }

    ## Assigning labels to the tick marks
    if (!is.null(xtlab)) 
    {
        xLabels <- xtlab
    } else {
        xLabels <- as.character(xaxisTicks)
    }

    ## Avoiding too many tick marks
    if (xttrim)
    {
        lenXT <- length(xaxisTicks)
        if (lenXT > 6) 
        {
            thinFactor <- max(c(2, floor(lenXT/6)))
            halfLXT <- floor(lenXT / thinFactor) - 1
            chosenInd <- 1 + thinFactor*(0:halfLXT)  
            # "1" is ensuring that control always is present
            xaxisTicks <- xaxisTicks[chosenInd]
            xLabels <- xLabels[chosenInd]
        }
    }     

    ## Assigning special name to first tick mark
    if (logX && (is.null(xtlab)) && (!is.null(conName)))
    {
        xLabels[1] <- conName
    }   

    ## Updating labels
    xLabels <- as.expression(xLabels)

    ## Updating x axis labels
#    xLabels <- as.expression(xaxisTicks)
#    if (conNameYes) {xLabels[1] <- conName}                                                

    axis(1, at = xaxisTicks, labels = xLabels, cex.axis = cex.axis)
    axis(2, at = yaxisTicks, labels = yLabels, cex.axis = cex.axis)        
}



## Creating a broken axis
"brokenAxis" <- function(bcontrol, broken, bp, dosePts, gridsize, log, logX, logDose, plotit = TRUE)
{
    notNullLD <- !is.null(logDose) 
    if ((broken) && (logX  || (notNullLD)))
    {
        bList <- list(factor = 2, style = "slash", width = 0.02)
            
        if (!is.null(bcontrol))
        {
            namesBC <- names(bcontrol)
            for (j in 1:length(bcontrol))
            {
                bList[[namesBC[j]]] <- bcontrol[[j]]
            }
        }
        breakStyle <- bList$"style"  # "slash"
        breakWidth <- bList$"width"  # 0.02  # default in axis.break
        clFactor <- bList$"factor"  # 2
         
        if (notNullLD)
        {
            brokenx <- log(clFactor * (logDose^bp), logDose)
        } else {
            brokenx <- clFactor * bp
        }             
#        brokenx <- clFactor * bp                   
        if ( (log == "x") || (log == "xy") || (log == "yx") )
        {
            ivMid <- dosePts > brokenx            
        } else {
            ivMid <- rep(TRUE, gridsize)
        }
        if (plotit)
        {
#            require(plotrix, quietly = TRUE)
            axis.break(1, brokenx, style = breakStyle, brw = breakWidth)        
        }
        
    } else {
        ivMid <- rep(TRUE, gridsize)
    }
    return(ivMid)
}


## Adding legend and legend text
"makeLegend" <- function(colourVec, legend, legendCex, legendPos, legendText, lenlev, level, ltyVec, noPlot, pchVec, type,
xLimits, yLimits)
{
    if (!legend) {return(invisible(NULL))}
    
    legendLevels <- as.character(level)
    if (!missing(legendText)) 
    {
        lenLT <- length(legendText)
    
        if (lenLT == lenlev) {legendLevels <- legendText}
        
        if (lenLT == 1) {legendLevels <- rep(legendText, lenlev)}
    }
    levInd <- 1:lenlev
    
    ## Removing line types when lines are not drawn
    ltyVec[noPlot] <- 0    
    if (identical(type, "obs"))
    {
        ltyVec[levInd] <- 0
    }
    
    ## Removing plot symbol when no points are drawn
    if ( (identical(type, "none")) || (identical(type, "bars")) )
    {
        pchVec[levInd] <- NA
    }
      
    ## Defining position of legend
    if (!missing(legendPos))
    {
        if ( (is.numeric(legendPos)) && (length(legendPos) == 2) )
        xlPos <- legendPos[1]
        ylPos <- legendPos[2]
    } else {
        xlPos <- xLimits[2] 
        ylPos <- yLimits[2]
    }
    
    ## Adding legend
    legend(xlPos, ylPos, legendLevels, lty = ltyVec[levInd], pch = pchVec[levInd], 
    col = colourVec[levInd], bty = "n", xjust = 1, yjust = 1, cex = legendCex)
}

"parFct" <- function(gpar, lenlev, defVal = NULL)
{
    if (!missing(gpar)) 
    {
        if (length(gpar) == 1) 
        {
            return(rep(gpar, lenlev))
        } else {
            return(gpar)
        }
    } else {
        if (is.null(defVal)) {return(1:lenlev)} else {rep(defVal, lenlev)}
    }
}


getLU <- function(object)
{
    parmMat <- object$"parmMat"
#    rownames(parmMat) <- object$"parNames"[[2]]

    fixedVal <- object$fct$fixed
    lenFV <- length(fixedVal)
#    parmMatExt <- matrix(NA, length(fixedVal), ncol(parmMat))
    parmMatExt <- matrix(fixedVal, length(fixedVal), ncol(parmMat))
    parmMatExt[is.na(fixedVal), ] <- parmMat

    parmMatExt
}


normalizeLU <- function(x, y, normRef = 1)
{
    cVal <- y[2]
    dVal <- y[3]  
    normRef * ((x - cVal) / (dVal - cVal))
}

#mapply(normalizeLU, as.list(as.data.frame(matrix(1:20, 10, 2))), as.list(as.data.frame(getLU(S.alba.m1))))

