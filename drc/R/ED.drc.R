ED <- function (object, ...) UseMethod("ED", object)

"ED.drc" <-
function(object, respLev, interval = c("none", "delta", "fls", "tfls"), clevel = NULL,
level = ifelse(!(interval == "none"), 0.95, NULL), reference = c("control", "upper"), 
type = c("relative", "absolute"), lref, uref, bound = TRUE, od = FALSE, vcov. = vcov, # robust = false,
display = TRUE, pool = TRUE, logBase = NULL, multcomp = FALSE, ...)
{
    interval <- match.arg(interval)
    reference <- match.arg(reference)
    type <- match.arg(type)
    
    ## Checking 'respLev' vector ... should be numbers between 0 and 100
    if ( (type == "relative") && (bound) ) 
    {
        if (any(respLev <= 0 | respLev >= 100)) 
        {
            stop("Response levels (percentages) outside the interval ]0, 100[ not allowed")
        }
    }

    ## Retrieving relevant quantities
    EDlist <- object$fct$"edfct"  
    if (is.null(EDlist)) {stop("ED values cannot be calculated")}         
    indexMat <- object$"indexMat"
#    parmMat <- matrix(coef(object)[indexMat], ncol = ncol(indexMat))   
    parmMat <- object$"parmMat"
#    strParm0 <- colnames(parmMat)    
    strParm0 <- sort(colnames(object$"parmMat"))
    
    curveNames <- colnames(object$"parmMat")
    options(warn = -1)  # switching off warnings caused by coercion in the if statement
    if (any(is.na(as.numeric(curveNames))))
    {
        curveOrder <- order(curveNames)
    } else { # if names are numbers then skip re-ordering
        curveOrder <- 1:length(curveNames)
    }
    options(warn = 0)  # normalizing behaviour of warnings
    
    strParm0 <- curveNames[curveOrder]
    indexMat <- indexMat[, curveOrder, drop = FALSE]
    parmMat <- parmMat[, curveOrder, drop = FALSE]
    
    strParm <- strParm0
#    if (robust)
#    {
#        vcMat <- sandwich(object)
#    } else {
#        vcMat <- vcov(object, od = od, pool = pool)
#    }
    vcMat <- vcov.(object)
    
    ## Defining vectors and matrices
    ncolIM <- ncol(indexMat)
    indexVec <- 1:ncolIM    
#    lenEB <- ncolIM    
    lenPV <- length(respLev)  # used twice below
    noRows <- ncolIM * lenPV
    dimNames <- rep("", noRows)  # lenEB*lenPV, 2)
    EDmat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)
    oriMat <- matrix(0, noRows, 2)  # lenEB*lenPV, 2)

    ## Skipping curve id if only one curve is present
#    lenIV <- lenEB  # ncol(indexMat)
    if (identical(length(unique(strParm)), 1)) 
    {
#        strParm[1:lenIV] <- rep("", lenIV)
        strParm[indexVec] <- rep("", ncolIM)
    } else {
        strParm <- paste(strParm, ":", sep = "")
    }

    ## Calculating estimates and estimated standard errors
    rowIndex <- 1
    lenIV <- length(indexVec)
    dEDmat <- matrix(0, lenPV * lenIV, nrow(vcMat))
    for (i in indexVec)
    {
        parmChosen <- parmMat[, i]
#        print(parmChosen)
        parmInd <- indexMat[, i]
        varCov <- vcMat[parmInd, parmInd]

        if ((is.null(clevel)) || (strParm0[i] %in% clevel))
        {
        for (j in 1:lenPV)
        {
            EDeval <- EDlist(parmChosen, respLev[j], reference = reference, type = type, ...)            
            EDval <- EDeval[[1]]
            dEDval <- EDeval[[2]]
#            print(c(i,j,parmInd))
#            print(dEDval)
            dEDmat[(i-1)*lenPV + j, parmInd] <- dEDval 
            
            oriMat[rowIndex, 1] <- EDval
            oriMat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)
                   
            if (!is.null(logBase))
            {
                EDval <- logBase^(EDval)                
                dEDval <- EDval * log(logBase) * dEDval
            }
            EDmat[rowIndex, 1] <- EDval
            EDmat[rowIndex, 2] <- sqrt(dEDval %*% varCov %*% dEDval)

            dimNames[rowIndex] <- paste(strParm[i], respLev[j], sep = "")
            rowIndex <- rowIndex + 1
        }
        } else {
            rowsToRemove <- rowIndex:(rowIndex + lenPV - 1)
            EDmat <- EDmat[-rowsToRemove, , drop = FALSE]
            dimNames <- dimNames[-rowsToRemove]
        }
    }
    
    ## Defining column names
    colNames <- c("Estimate", "Std. Error")
    
#    ## Using t-distribution for continuous data
#    ## (only under the normality assumption)
#    if (object$"type" == "continuous")
#    {
#        qFct <- function(x) {qt(x, df.residual(object))}
#    } else { # Otherwise the standard normal distribution is used
#        qFct <- qnorm
#    }

    ## Calculating the confidence intervals
    if (interval == "delta")
    {
#        colNames <- c(colNames, "Lower", "Upper")    
#        ciMat <- matrix(0, lenEB*lenPV, 2)
#        tquan <- qFct(1 - (1 - level)/2)        
#        ciMat[, 1] <- EDmat[, 1] - tquan * EDmat[, 2]
#        ciMat[, 2] <- EDmat[, 1] + tquan * EDmat[, 2]
        
#        print(EDmat)
        intMat <- confint.basic(EDmat, level, object$"type", df.residual(object), FALSE)
#        print(intMat)
        intLabel <- "Delta method"
    }
    
    if (interval == "tfls")
    {
#        colNames <- c( colNames, "Lower", "Upper")     
#        lsVal <- log(oriMat[, 1])
#        lsdVal <- oriMat[, 2]/oriMat[, 1]
#        ciMat <- matrix(0, lenEB*lenPV, 2)
#        tquan <- qFct(1 - (1 - level)/2)                       
#        ciMat[, 1] <- exp(lsVal - tquan * lsdVal)
#        ciMat[, 2] <- exp(lsVal + tquan * lsdVal)
        
        intMat <- exp(confint.basic(matrix(c(log(oriMat[, 1]), oriMat[, 2] / oriMat[, 1]), ncol = 2), 
        level, object$"type", df.residual(object), FALSE))
        intLabel <- "To and from log scale"       
    }

    if (interval == "fls")
    {        
#        ciMat <- matrix(0, lenEB*lenPV, 2)
#        tquan <- qFct(1 - (1 - level)/2) 
        if (is.null(logBase)) 
        {
            logBase <- exp(1)
            EDmat[, 1] <- exp(EDmat[, 1])  # back-transforming log ED values
        }
#        ciMat[, 1] <- logBase^(oriMat[, 1] - tquan * oriMat[, 2])
#        ciMat[, 2] <- logBase^(oriMat[, 1] + tquan * oriMat[, 2])

        intMat <- logBase^(confint.basic(oriMat, level, object$"type", df.residual(object), FALSE))
        intLabel <- "Back-transformed from log scale"  

        ## Dropping estimated standard errors (not relevant after back transformation)        
        EDmat <- EDmat[, -2, drop = FALSE]  
        colNames <- colNames[-2]
#        colNames <- c(colNames[-2], "Lower", "Upper")  # standard errors not relevant        
    }
    
    if (identical(interval, "none"))
    {
        intLabel <- NULL
    } else {
        EDmat <- as.matrix(cbind(EDmat, intMat))
        colNames <- c(colNames, "Lower", "Upper")         
    }   
    dimnames(EDmat) <- list(dimNames, colNames)
    resPrint(EDmat, "Estimated effective doses", interval, intLabel, display = display)
    
    ## require(multcomp, quietly = TRUE)
#    invisible(list(EDdisplay = EDmat, EDmultcomp = list(EDmat[, 1], dEDmat %*% vcMat %*% t(dEDmat))))  
    
    if(multcomp)
    {  
        invisible(list(EDdisplay = EDmat, 
                       EDmultcomp = parm(EDmat[, 1], (dEDmat %*% vcMat %*% t(dEDmat))[1:nrow(EDmat), 1:nrow(EDmat), drop = FALSE])))
    } else {
        invisible(EDmat)
      
    }   
}

