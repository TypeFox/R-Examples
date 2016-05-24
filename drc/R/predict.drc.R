"predict.drc" <- function(object, newdata, se.fit = FALSE, interval = c("none", "confidence", "prediction"), 
level = 0.95, na.action = na.pass, od = FALSE, vcov. = vcov, ...)
{
    ## Checking arguments
    interval <- match.arg(interval)
    respType <- object[["type"]]

    dataList <- object[["dataList"]]    
    doseDim <- ncol(dataList[["dose"]])
    if (is.null(doseDim)) {doseDim <- 1} 

    ## Assigning dataset from object if no data frame is provided
    if (missing(newdata)) 
    {
#        predValues <- fitted(object)  # not used
#        newdata <- data.frame(object$data[, 1], object$data[, 3])
#        dataList <- object[["dataList"]]

        ## New part (25/6-2014)
        doseVec <- dataList[["dose"]]
        if (identical(respType, "event"))
        {
            groupLevels <- as.character(dataList[["plotid"]])
        } else {
            groupLevels <- as.character(dataList[["curveid"]])
        }        
#                
#        if (identical(respType, "event"))
#        {
#            newdata <- data.frame(dataList[["dose"]], dataList[["plotid"]])
#        } else {
#            newdata <- data.frame(dataList[["dose"]], dataList[["curveid"]])
#        }
    } else {
        if (ncol(newdata) < (doseDim + 1)) {newdata <- data.frame(newdata, rep(1, nrow(newdata)))}
#        ndncol <- ncol(newdata)
#        doseVec <- newdata[, 1:(ndncol-1)]
        doseVec <- newdata[, 1:doseDim]
#        groupLevels <- as.character(newdata[, ndncol])  # 'as.character()' used to suppress factor levels         
        groupLevels <- as.character(newdata[, doseDim + 1])  # 'as.character()' used to suppress factor levels         
    }
    noNewData <- length(groupLevels)
    
#    if (ncol(newdata) < 2) {newdata <- data.frame(newdata, rep(1, nrow(newdata)))}
#    if (ncol(newdata) > 2) {stop("More than 2 variables in 'newdata' argument")}
    
    ## Defining dose values -- dose in the first column!   
#    doseVec <- newdata[, 1]
#    groupLevels <- as.character(newdata[, 2])  # 'as.character()' used to suppress factor levels 
#    noNewData <- length(doseVec)


    ## Transforming to dose scale if necessary
    powerExp <- (object$"curve")[[2]]
    if (!is.null(powerExp))
    {
        doseVec <- powerExp ^ doseVec
    }

    ## Retrieving matrix of parameter estimates
    parmMat <- object[["parmMat"]] 
    pm <- t(parmMat[, groupLevels, drop = FALSE])
           
#    parmNames <- colnames(parmMat)
#    lenCN <- length(parmNames)
#    indVec <- 1:lenCN
#    names(indVec) <- parmNames
#    if (lenCN > 1)
#    {
#        indVec <- indVec[as.character(newdata[, 2])]
#        
##        groupLevels <- newdata[, 2]
#        if (!all(is.numeric(groupLevels)))
#        {
##            pm <- parmMat[, as.character(groupLevels)]  # 'as.character()' used to suppress factor levels            
#            pm <- parmMat[, groupLevels]            
#        } else {
#            pm <- parmMat[, groupLevels]
#        }
#        pm <- parmMat[, groupLevels]
#        
#    } else {
#        lenDV <- length(doseVec)
##        indVec <- rep(1, lenDV)
#        pm <- matrix(parmMat[, 1], length(parmMat[, 1]), lenDV)
#    }    

#    ## Checking for NAs in matrix of parameter estimates
#    naVec <- rep(FALSE, lenCN)
#    for (i in 1:lenCN)
#    {
#        naVec[i] <- any(is.na(parmMat[, i]))
#    }
#    parmMat <- parmMat[, !naVec, drop = FALSE] 


    ## Retrieving variance-covariance matrix
    sumObj <- summary(object, od = od)
#    varMat <- sumObj[["varMat"]]  
    vcovMat <- vcov.(object)      

    ## Defining index matrix for parameter estimates
    indexMat <- object[["indexMat"]]
    
    ## Calculating predicted values  
#    indexVec <- as.vector(indVec)  
#    print(indexVec)  
#    lenIV <- length(indexVec)    
    
    
#    retMat <- matrix(0, lenIV, 4)
    retMat <- matrix(0, noNewData, 4)
    colnames(retMat) <- c("Prediction", "SE", "Lower", "Upper")
    objFct <- object[["fct"]]
    retMat[, 1] <- objFct$"fct"(doseVec, pm)
#    print(pm)
    
    ## Checking if derivatives are available
    deriv1 <- objFct$"deriv1"
    if (is.null(deriv1))
    {
        return(retMat[, 1])        
    }    

    ## Calculating the quantile to be used in the confidence intervals
    if (!identical(interval, "none"))
    {    
        if (identical(respType, "continuous"))
        {
            tquan <- qt(1 - (1 - level)/2, df.residual(object))   
        } else {
            tquan <- qnorm(1 - (1 - level)/2)
        }
    }  
    
    ## Calculating standard errors and/or confidence intervals
    if (se.fit || (!identical(interval, "none")))
    {
         if (identical(interval, "prediction"))
         {
             sumObjRV <- sumObj$"resVar"
         } else {
             sumObjRV <- 0
         }
#        rowIndex <- 1    
#        for (i in indexVec)
#        for (i in 1:ncol(indexMat))

#        groupLevels <- newdata[, 2]
        piMat <- indexMat[, groupLevels, drop = FALSE]
#        print(piMat)
#        print(groupLevels)
        for (rowIndex in 1:noNewData)
        {
#            parmInd <- indexMat[, i]
#            print(indexVec) 
#            print(varMat)
#            print(parmInd)       
            
#            varCov <- varMat[parmInd, parmInd]
#            print(varCov)
#            groupLevels <- newdata[, 2]
#            parmInd <- indexMat[, groupLevels[rowIndex]]
#            varCov <- varMat[parmInd, parmInd]

            parmInd <- piMat[, rowIndex] 
            varCov <- vcovMat[parmInd, parmInd]

#            parmChosen <- t(parmMat[, i, drop = FALSE])
#            parmChosen <- t(pm[, rowIndex, drop = FALSE])
#            dfEval <- deriv1(doseVec[rowIndex], parmChosen)

            dfEval <- deriv1(doseVec[rowIndex], pm[rowIndex, , drop = FALSE])
            varVal <- dfEval %*% varCov %*% dfEval
            retMat[rowIndex, 2] <- sqrt(varVal)  
#            retMat[rowIndex, 2] <- sqrt(dfEval %*% varCov %*% dfEval)  

            if (!se.fit)
            {
                retMat[rowIndex, 3:4] <- retMat[rowIndex, 1] + (tquan * sqrt(varVal + sumObjRV)) * c(-1, 1)
#                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal + sumObjRV)
#                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal + sumObjRV)   
            }    
#            if (identical(interval, "confidence"))
#            {
#                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal)
#                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal)            
#            }
#            if (identical(interval, "prediction"))
#            {
#                sumObjRV <- sumObj$"resVar"
#                retMat[rowIndex, 3] <- retMat[rowIndex, 1] - tquan * sqrt(varVal + sumObjRV)
#                retMat[rowIndex, 4] <- retMat[rowIndex, 1] + tquan * sqrt(varVal + sumObjRV)                        
#            }          
#            rowIndex <- rowIndex + 1        
        }
    }
    ## Keeping relevant indices
    keepInd <- 1
    if (se.fit) {keepInd <- c(keepInd, 2)}
    if (!identical(interval, "none")) {keepInd <- c(keepInd, 3, 4)}
    
    return(retMat[, keepInd])  # , drop = FALSE])
}


