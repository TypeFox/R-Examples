"relpot" <- function(object, plotit = TRUE, compMatch = NULL, percVec = NULL, interval = "none", 
type = c("relative", "absolute"), scale = c("original", "percent", "unconstrained"), ...)
{
    scale <- match.arg(scale)
    type <- match.arg(type)

#    ## Checking arguments
#    if (length(compMatch) != 2)
#    {
#        stop("Argument 'compMatch' should have length 2")
#    }

    ## Defining range for 'percVec' 
    parmMat <- commatFct(object, compMatch)    
    lowerVec <- apply(parmMat, 2, object$"fct"$"lowerAs")
    upperVec <- apply(parmMat, 2, object$"fct"$"upperAs")
    maxLow <- max(lowerVec)
    minUp <- min(upperVec)
    
    if ( (type == "absolute") && (is.null(percVec)) )
    {
        percVec <- seq(maxLow*1.05, minUp*0.95, length.out = 100)      
    }
    if ( (type == "relative") && (is.null(percVec)) )
    {
        uMin <- max( (maxLow - lowerVec) / (upperVec - lowerVec) )
        uMax <- min( (minUp - lowerVec) / (upperVec - lowerVec) )  
        if (object$"fct"$"monoton"(parmMat[, 1]) < 0)
        {
            uTemp <- uMin
            uMin <- 1 - uMax
            uMax <- 1 - uTemp
        }  
        percVec <- 100 * (seq(uMin, uMax, length.out = 101))[-c(1, 101)]
    }
    if ( (type == "relative") && (scale == "unconstrained") )
    {
        percVec <- 1:99
    } 

    lenpv <- length(percVec)
    rpVec <- rep(NA, lenpv)
    if (identical(interval, "none"))
    {
        for (i in 1:lenpv)
        {
            SIobj <- SI(object, rep(percVec[i], 2), compMatch, type = type, display = FALSE)[["SIdisplay"]]
            rpVec[i] <- SIobj[1]
        }
    } else {
        lrpVec <- rep(NA, lenpv)
        urpVec <- rep(NA, lenpv)

        for (i in 1:lenpv)
        {
            SIobj <- SI(object, rep(percVec[i], 2), compMatch, interval = interval, type = type, display = FALSE)[["SIdisplay"]]
            rpVec[i] <- SIobj[1]
            lrpVec[i] <- SIobj[2]
            urpVec[i] <- SIobj[3]            
        }
    }

    if (plotit)
    {
        if ( (type == "relative") && ((scale == "percent") || (scale == "unconstrained")) )
        {
            xlabStr <- "Relative response level (%)"
            xVec <- percVec
        }
        if ( (type == "relative") && (scale == "original") )
        {
            xlabStr <- "Response level"
            xVec <- seq(maxLow, minUp, length.out = 99)            
        }
        if (type == "absolute")
        {
            xlabStr <- "Response level"
            xVec <- percVec
        }
                
        if (!identical(interval, "none"))
        {
            plot(xVec, rpVec, type = "l", xlab = xlabStr, ylab = "Relative potency", 
            ylim = c(min(lrpVec), max(urpVec)), ...)
         
            lines(xVec, lrpVec, lty = 3)
            lines(xVec, urpVec, lty = 3)            
        } else {
            plot(xVec, rpVec, type = "l", xlab = xlabStr, ylab = "Relative potency", ...)
        }

        ## Adding reference line corresponding to EC50/ED50
        if (type == "relative")
        { 
            abline(h = SI(object, c(50, 50), compMatch, type = type, display = FALSE)[["SIdisplay"]][1], lty = 2)
        }
        
    }
    invisible(list(x = xVec, y = rpVec, percVec = percVec))
}
