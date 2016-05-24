"estfun.drc" <- function (x, ...)
#"estfun.drc" <- function (x, cvar = NULL, ...)  
{
    ## Extending the matrix of derivatives to have one column per parameter (not only per model parameter)
    ##  only relevant in case several curves (with different parameters) were fitted
    xderiv1 <- x$deriv1
    indexMat0 <- x$indexMat2
    if (is.matrix(indexMat0))
    {
        indMat <- t(indexMat0)
    } else {
        indMat <- t(matrix(indexMat0, nrow = 1))
    }
#    print(indMat)
    colnames(indMat) <- colnames(x$indexMat)
    curveID <- x$dataList[["curveid"]]
#    xderiv2 <- xderiv1[, rep(1:ncol(xderiv1), apply(indMat, 1, length))]

#    xderiv2 <- xderiv1[, rep(1:ncol(xderiv1), apply(indMat, 1, function(x){length(unique(x))}))]
#    cnInd <- colnames(indMat)
#    for (i in 1:ncol(indMat))
#    {
#        xderiv2[curveID != cnInd[i], indMat[, i]] <- 0
#    }
##    xderiv2[curveID != colnames(indMat)[2], indMat[, 2]] <- 0
##    print(xderiv2)

    ## Defining a helper function for making a wide version of the mattrix of parameter estimates
    xderiv2Fct <- function(xderiv1, indMat, curveID)
    {
        xderiv2 <- xderiv1[, rep(1:ncol(xderiv1), apply(indMat, 1, function(x){length(unique(x))}))]    
#        print(xderiv2)
        
        cnInd <- colnames(indMat)
        for (i in 1:ncol(indMat))
        {
            xderiv2[curveID != cnInd[i], indMat[, i]] <- 0
        }
        xderiv2
    }

    if (identical(x$type, "continuous"))
    {
        xderiv2 <- xderiv2Fct(xderiv1, indMat, curveID)    
#        rval <- (weights(x) * residuals(x)) * x$deriv1
        rval <- (weights(x) * residuals(x)) * xderiv2
    }
    
    if (identical(x$type, "binomial"))
    {  # not handling missing values
        nTotal <- weights(x)
        nObs <- x[["dataList"]][["resp"]] * nTotal
        fittedVal <- fitted(x)
        rval0 <- nObs/fittedVal + (nObs - nTotal)/(1 - fittedVal)
#        rval <- x$deriv1 * rval0
        xderiv2 <- xderiv2Fct(xderiv1, indMat, curveID)
        rval <- xderiv2 * rval0    
        rval
    }
    if (identical(x$type, "Poisson"))
    {  # not handling missing values
        resp <- x[["dataList"]][["resp"]]
        fittedVal <- fitted(x)
        rval0 <- resp / fittedVal - 1
#        rval <- x$deriv1 * rval0
        xderiv2 <- xderiv2Fct(xderiv1, indMat, curveID)  
        rval <- xderiv2 * rval0  
        rval
    }    
    if (identical(x$type, "event"))
    {  # not handling missing values
    
        Ft1 <- predict(x, data.frame(x$data[, 1], x$data[, 5]))
        Ft2 <- predict(x, data.frame(x$data[, 2], x$data[, 5]))
        Ft2[!is.finite(x$data[, 2])] <- 1    
        diffF <- Ft2 - Ft1
    
        DFt1.0 <- x[["fct"]]$"deriv1"(x$data[, 1], t(x[["parmMat"]][, as.character(x$data[, 5])]))        
        DFt2.0 <- x[["fct"]]$"deriv1"(x$data[, 2], t(x[["parmMat"]][, as.character(x$data[, 5])]))
        DFt2.0[!is.finite(x$data[, 2]), ] <- 0         
        diffDF.0 <- DFt2.0 - DFt1.0 
        diffDF <- xderiv2Fct(diffDF.0, indMat, x$data[, 5])
         
        rval <- (x$data[, 3] / diffF) * diffDF  
        
        if (FALSE)
        { 
        dataList <- x[["data"]]
        xderiv1 <- x[["fct"]]$"deriv1"(dataList[, 1], t(x[["parmMat"]][, as.character(dataList[, 4])]))
        print(xderiv1)
        xderiv2 <- xderiv2Fct(xderiv1, indMat, data[, 4])
#        print(xderiv2)
        resp <- c(0, dataList[, 3])
        fittedVal <- c(0, predict(dataList[, 2]))
#        print(resp)
#        print(fittedVal)
#        fittedVal[length(fittedVal)] <- 1
#        fittedVal2 <- c(fittedVal[-1], 1)  # assuming data ordered according to time
#        rval0 <- diff(resp) / diff(fittedVal)
#        print(rval0)
##        lagDeriv1 <- apply(x$deriv1, 2, function(x){diff(x)})
#        lagDeriv1 <- apply(xderiv2, 2, function(x){diff(x)})
#        rval <- lagDeriv1 * rval0
#        rval <- xderiv2 * rval0    
#        rval
        rval <- xderiv2 * (diff(resp) / diff(fittedVal))
        }
    }
#    ## Summing up according to specified clusters (the variable "cvar")
#    if (!is.null(cvar))
#    {
#        rval <- rowsum(rval, cvar)
#    }    
        
#    print(rval)
    colnames(rval) <- names(coef(x))
#    colnames(rval) <- x$fct$names       
    rval
}

"bread.drc" <- function (x, ...) 
{
#    if (identical(x$type, "binomial"))
#    {
#        breadMat <- vcov(x) * unlist(x$sumList[1])
#    }
    if (identical(x$type, "continuous"))
    {
#        breadMat <- summary(x)$cov.unscaled * unlist(x$sumList[1])
        breadMat <- vcov(x) / (summary(x)$rse[1]^2) * unlist(x$sumList[1])
    } else { ## Note: not checked for event time data!
        breadMat <- vcov(x) * unlist(x$sumList[1])
    }
    return(breadMat)
}

