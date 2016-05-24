### Model-robust dose-response modelling

## Calculating leave-one-out predictions for the parametric and non-parametric fits separately
leaveOneOut <- function(object1, object2, dose, dataSet, resp, fixedEnd)
{   
    ## Leave-one-out predictions
    uniDose <- sort(unique(dose))
    lenUd <- length(uniDose)    
    
    doseDF <- dataSet[, tail(as.character(formula(object1)[[3]]), 1), drop = FALSE]  # picking dose column
    pred1 <- list()
    pred2 <- list()    
    for (i in 1:lenUd)
    {
        reFit1 <- update(object1, data = subset(dataSet, dose != uniDose[i]))
#        pred1[[i]] <- as.vector(predict(reFit1, newdata = subset(dataSet[, 1, drop = FALSE], dose == uniDose[i]), 
        pred1[[i]] <- as.vector(predict(reFit1, newdata = subset(doseDF, dose == uniDose[i]), 
        se.fit = FALSE))
#        print(pred1[[i]])
        
        reFit2 <- update(object2, data = subset(dataSet, dose != uniDose[i]))
#        control = loess.control(surface = "direct"))
        pred2[[i]] <- predict(reFit2, newdata = subset(doseDF, dose == uniDose[i]))
#        pred2[[i]] <- predict(reFit2, newdata = subset(dataSet[, 1, drop = FALSE], dose == uniDose[i]))
    }
    
    # Avoiding overflow problems
    pred2Vec <- as.vector(unlist(pred2))
    pred2Vec[pred2Vec < 0.01] <- 0.01  
    pred2Vec[pred2Vec > 0.99] <- 0.99
    
    ## Fixing boundary values at observed averages
    if (fixedEnd)
    {
        pred2Vec[1] <- mean(resp[dose == uniDose[1]])
        pred2Vec[lenUd] <- mean(resp[dose == uniDose[lenUd]])        
    }
    
    return(list(pred1 = as.vector(unlist(pred1)), pred2 = pred2Vec))
}


## Calculating weights to be used in PRESS* criterion under least squares
pressWeights <- function(w, lenData, nVec, object1, resp, object2)
{
    switch(w,
    "ad hoc" =  # similar to what Nottingham and Birch (2000) did
    {
        rVec <- resp * nVec
        any01 <- abs(rVec - nVec) < 1
        if (any(any01)) {rVec[any01] <- rVec[any01] - 0.5}
        pVec <- rVec / nVec
        nVec / (pVec * (1 - pVec))
    },
#    "inverse" = 
#    {
#        predVec1 <- as.vector(predict(object1, se.fit = FALSE))  # se.fit = FALSE not needed 
#        (predVec1 * (1 - predVec1)) / nVec
#    },    
    "none" = rep(1, lenData),            
    "nonpar" = 
    {
        predVec2 <- predict(object2)
        nVec / (predVec2 * (1 - predVec2))
    },      
    "par" = 
    {
        predVec1 <- as.vector(predict(object1, se.fit = FALSE))  # se.fit = FALSE not needed 
        nVec / (predVec1 * (1 - predVec1))
    },    
    "response" = nVec / (resp * (1 - resp)))
}

## Function calculating degrees of freedom
dfFct <- function(object1, object2)  # , trace1 = traceHat.drc, trace2 = traceHat.loess)
{
    # Trace hat for drm() model fit
    lenData <- object1$"sumList"$lenData
    traceHat.drc <- function(object)
    {
        lenData - df.residual(object)  # number of parameters ... for sure an easier way?
    }

    # Trace hat for loess fit
    traceHat.loess <- function(object)
    {
        object$trace.hat
    }
    
    function(lambda)
    {
#        lenUd - ( (1 - lambda) * traceHat.drc(object1) + lambda * traceHat.loess(object2))
         lenData - ( (1 - lambda) * traceHat.drc(object1) + lambda * traceHat.loess(object2))
    }
}

## Calculating leave-one-out predictions for the semi-parametric fit
predFct <- function(looList)
{
    function(lambda)
    {
        (1 - lambda) * looList$"pred1" + lambda * looList$"pred2"
    }
}


## Obtaining model-robust fit (evaluating PRESS* criterion)
"mrdrm" <- function(object1, object2, lambda = (0:10)/10, criterion = c("gcv", "lcv"), critFct = c("ls", "ll"),
ls.weights = c("nonpar", "ad hoc", "none", "par", "response"), fixedEnd = FALSE, unitScale = FALSE)  
# object1 is the drm() model fit
# fixedEnd = TRUE favours the non-parametric model!!!
{
    criterion <- match.arg(criterion)
    critFct <- match.arg(critFct)    
    ls.weights <- match.arg(ls.weights)
    ## Enforcing the least squares criterion function without weighting for continuous data
    if (identical(object1$"type", "continuous")) 
    {
        critFct <- "ls"
        ls.weights <- "none"
    } 

    ## Fitting a local linear regression with default settings in case no fit is provided
    if (missing(object2))
    {
        object2 <- loess(formula(object1), data = object1$"origData", degree = 1)
    }

    ## Checking that a local linear regression fit is supplied
    if ( (object2$"pars"$"degree" > 1) || (object2$"pars"$"degree" < 1) ) 
    {
        stop("Local regression fit not linear!", call. = FALSE)
    }     

    ## Retrieving data
    dataSet <- object1$"origData"
    dataSet2 <- object1$"data"
    dose <- dataSet2[, 1]
#    uniDose <- unique(dose)
#    lenUd <- length(uniDose)    
#    lenData <- length(dose)  
    resp <- dataSet2[, 2]    
    lenData <- object1$"sumList"$lenData

    ## Transforming doses into the unit interval
    if (unitScale)
    {
        uniqDose <- sort(unique(dose))
        lenUD <- length(uniqDose)
        doseLoess <- loess((0:(lenUD - 1))/lenUD ~ uniqDose)
        dosePredict <- function(dose) {predict(doseLoess, data.frame(uniqDose = dose))}
        unitDose <- dosePredict(dose)

        object2 <- loess(resp ~ unitDose, degree = 1)          
    } else {
        dosePredict <- function(dose) {dose}  # identity map
    }

       
#    predVec <- lambda * looList$"pred1" + (1 - lambda) * looList$"pred2"
##    predVec <- lambda * as.vector(unlist(pred1)) + (1 - lambda) * as.vector(unlist(pred2))
##    print(predVec)
    
    ## Press value
    dFct <- dfFct(object1, object2)    
    nVec <- object1$weights
    pressFct <- switch(critFct,
    "ls" =  # least squares criterion function
    {
        ## Weights
#        varVec <- weightFct(w, lenData, nVec, object1, resp) 
#        varVec <- switch(w,
#        "ad hoc" =  # similar to what Nottingham and Birch (2000) did
#        {
#            rVec <- resp * nVec
#            any01 <- abs(rVec - nVec) < 1
#            if (any(any01)) {rVec[any01] <- rVec[any01] - 0.5}
#            pVec <- rVec / nVec
#            nVec / (pVec * (1 - pVec))
#        }, 
#        "model-based" = 
#        {
#            predVec1 <- as.vector(predict(object1, se.fit = FALSE))  
#            nVec / (predVec1 * (1 - predVec1))
#        }, 
#        "none" = rep(1, lenData),   
#        "response" = nVec / (resp * (1 - resp)))
##        wVec <- object1$"weights" / (resp * (1 - resp))
##        wVec <- object1$"weights" / varVec
#        print(varVec)

        ## Degrees of freedom
#        dfVal <- lenUd - (lambda * traceHat.drc(object1) + (1 - lambda) * traceHat.loess(object2))
    
        ## Press value
#        pressFct1(resp - predVec, varVec, dfVal)  
        
        switch(criterion,
        "gcv" = {  ## Using GCV
            looList <- NULL
            pFct <- predFct(list(pred1 = predict(object1), pred2 = predict(object2)))
            pwVec <- NULL  # No weights used
                
            function(lambda)
            {
                lenData * sum((resp - pFct(lambda))^2) / (dFct(lambda)^2)
            }
        },
        "lcv" = {  # Using CV
            ## Calculating leave-one-out predictions
            looList <- leaveOneOut(object1, object2, dose, dataSet, resp, fixedEnd)
            pFct <- predFct(looList)    
            pressFct1 <- function(r, w, den) {sum(w * (r^2) / den, na.rm = TRUE)} 
            pwVec <- pressWeights(ls.weights, lenData, nVec, object1, resp, object2)                
        
            function(lambda)
            {
                pressFct1(resp - pFct(lambda), pwVec, dFct(lambda))
            }
        })  
    },
    "ll" =  # log likelihood criterion function   
    {
        ## Calculating leave-one-out predictions
        looList <- leaveOneOut(object1, object2, dose, dataSet, resp, fixedEnd)
        pFct <- predFct(looList)
        pwVec <- NULL  # No weights used
        pressFct2 <- function(n, x, p) {-sum((n - x) * log(1 - p) + x * log(p), na.rm = TRUE)}  
        function(lambda)
        {
            pressFct2(nVec, resp * nVec, pFct(lambda))
        }
    })    
#    sum ( wVec * (resp - predVec)^2 / dfVal, na.rm = TRUE )  # lenUd - trace()

    vPressFct <- Vectorize(pressFct, "lambda")
    
    if (length(lambda) > 1)
    {
        pVec <- vPressFct(lambda)
        optimalLambda <- lambda[which.min(pVec)] 
    } else {
        pVec <- NA  # no PRESS criterion used
        optimalLambda <- lambda
    }
    
    ## Calculating fitted values
    pred1 <- as.vector(predict(object1))   
    pred2 <- predict(object2)      
    fitVal <- (1 - optimalLambda) * pred1 + optimalLambda * pred2
    
    ## Calculating goodness-of-fit values
#    gofVal <- switch(object1$"type",
#    "binomial" = c(sum ( nVec * ((resp - fitVal)^2) / (fitVal * (1 - fitVal))),
#    sum ( nVec * ((resp - pred1)^2) / (pred1 * (1 - pred1)))),  
#    # adjustment for fitVal close to 0 or 1?
#    "continuous" = c(sum( (resp - fitVal)^2 ), sum( (resp - pred1)^2 ))) 
#        
#    ## Calculating AIC values    
#    aicVal <- switch(object1$"type",
#    "binomial" = -sum((nVec - nVec * resp) * log(1 - fitVal) + (nVec * resp) * log(fitVal), na.rm = TRUE),
#    "continuous" = lenData * log(2*pi) + lenData * log(gofVal/lenData) + lenData + 2 * dFct(optimalLambda))    
#    
#    ## Calculating residuals standard error
#    seVal <- switch(object1$"type",
#    "binomial" = NA,
#    "continuous" = gofVal / dFct(optimalLambda))    
    
    dfVal <- dFct(optimalLambda)
    gofVec <- switch(object1$"type",
    "binomial" = {
        success <- nVec * resp
    
        c(sum(nVec * ((resp - fitVal)^2) / (fitVal * (1 - fitVal))),
        sum(nVec * ((resp - pred1)^2) / (pred1 * (1 - pred1))),
        -2 * sum(log(choose(nVec, success))) - 2 * sum((nVec - success) * log(1 - fitVal) + success * log(fitVal), 
        na.rm = TRUE) + 2 * (lenData - dfVal), NA)
    },
    "continuous" = {
        gofVal <- sum( (resp - fitVal)^2 )
        
        c(gofVal, sum( (resp - pred1)^2 ),
        lenData * log(2*pi) + lenData * log(gofVal/lenData) + lenData + 2 * (lenData - dfVal),
        gofVal / dfVal)
        # Cleveland (1979) generalised
    })
    names(gofVec) <- c("mr.gof", "p.gof", "aic", "rv")
        
    retList <- list(pressVal = pVec, lambda = optimalLambda, fitted = fitVal, gof = gofVec,
    object1 = object1, object2 = object2, dose = dose, EDmethod = "inverse", ll = looList,
    ls.weights = pwVec, df = dfVal)
    class(retList) <- "mrdrc"
    
    retList
}

## Calculating loess fit (not used, just to understand)
"loessEst" <- function(x0, x, y, span, logScale = FALSE)
{
    tricubic <- function(x, maxDist)
    {
        (1 - abs(x / maxDist)^3)^3
    }
    
        nsDistVec <- abs(x0 - x) 
        distVec <- sort(nsDistVec)
        
        x <- x[order(nsDistVec)]
        y <- y[order(nsDistVec)]        
       
        wVec <- rep(0, length(x))
        iVec <- distVec > quantile(distVec, span, type = 3)

        maxDist <- max(distVec[!iVec])  # max(distVec)
        if (!logScale)
        {
            wVec[!iVec] <- tricubic(distVec[!iVec], maxDist)
        } else {
            wVec[!iVec] <- tricubic(exp(distVec[!iVec]), exp(maxDist))
        }

    list(coef(lm(y ~ I(x - x0), weights = wVec))[1], wVec)
}

## Calculating hat matrix for loess fit
"hat.loess" <- function(x, span, x0 = x)
{
    tricubic <- function(x, maxDist)
    {
        (1 - abs(x / maxDist)^3)^3
    }
 
    X <- model.matrix( ~ x)
    X0 <- as.matrix(cbind(1, x0))
    lenx0 <- length(x0)
    lenx <- length(x)
    H <- matrix(0, lenx0, lenx)
    for (i in 1:lenx0)
    {
        distVec <- abs(x0[i] - x)
        wVec <- rep(0, lenx)
        iVec <- distVec < quantile(distVec, span)
        selectDistVec <- distVec[iVec]
        wVec[iVec] <- tricubic(selectDistVec, max(abs(selectDistVec)))
    
        Hi <- diag(wVec)
        H[i, ] <- X0[i, , drop = FALSE] %*% solve(t(X) %*% Hi %*% X) %*% t(X) %*% Hi
    }
    H
}

## Calculating hat matrix for drm() fit
"hat.drc" <- function(object, x, x0 = x)
{
    Dmat <- object$"deriv1"    

    pMat0 <- t(object$parmMat) 
    pMat <- matrix(pMat0, ncol = ncol(pMat0), nrow = length(x0), byrow = TRUE)  
    derivMat <- object$fct$deriv1(x0, pMat)
    
    switch(object$"type",
    "binomial" = {
        pred1 <- as.vector(predict(object, se.fit = FALSE))     
        
        # to avoid problems in 'Wmat' below
        pred1[pred1 < 0.01] <- 0.01  
        pred1[pred1 > 0.99] <- 0.99

        wMat <- diag(as.vector(object$"weights" / (pred1 * (1 - pred1))))
        derivMat %*% solve(t(Dmat) %*% wMat %*% Dmat) %*% t(Dmat) %*% wMat        
    },
    "continuous" =  derivMat %*% solve(t(Dmat) %*% Dmat) %*% t(Dmat))
}

## Calculating hat matrix for model-robust fit
"hat.mr" <- function(object, x0 = x)  # loess specific
{
    lambda <- object$"lambda"
    x <- object$"dose"
#    if (missing(x0)) {x0 <- x}
    
    (1 - lambda) * hat.drc(object$object1, x, x0) + lambda * hat.loess(x, object$object2$pars$span, x0)
}


"se.mr" <- function(object, x0)
{
    HMat <- hat.mr(object, x0)     

    switch(object$object1$"type",
    "binomial" = {

        pVec <- fitted(object)
        pVec[pVec < 0.01] <- 0.01  
        pVec[pVec > 0.99] <- 0.99
        ## Another option to use fitted values from parametric fit (for sure between 0 and 1)        

        varY <- diag(as.vector((pVec * (1 - pVec)) / object$object1$"weights"))    
        sqrt(diag(HMat %*% varY %*% t(HMat)))
    },
    "continuous" = {
#        # Cleveland (1979) generalised
#        dfVal <- object$object1$"sumList"$lenData - sum(diag(HMat))  
#        sigma2 <- object$"gof" / dfVal
        sqrt(diag(HMat %*% t(HMat)) * object$"gof"[4])  # inefficient to do complete matrix multiplication
    })
}
 

## predict method for model-robust fit
predict.mrdrc <- function(object, newdata, se.fit = FALSE, interval = c("none", "confidence", "prediction"), 
level = 0.95, pava = FALSE, ...)
{
    interval <- match.arg(interval)
     
    if (missing(newdata))
    {
        predVec <- fitted(object)
        seVec <- se.mr(object, object$"dose")
    } else {
        pred1 <- as.vector(predict(object$"object1", newdata = newdata, se.fit = FALSE))    
        pred2 <- predict(object$"object2", newdata = newdata)      
        lambda <- object$"lambda"
        predVec <- (1 - lambda) * pred1 + lambda * pred2
        seVec <- se.mr(object, newdata[, 1])
    }
        
    if (se.fit)  # overrules the argument "interval"
    {
        retMat <- as.matrix(cbind(predVec, seVec))
        colnames(retMat) <- c("Prediction", "SE")
        return(retMat)
    }
    
    if (interval == "confidence")
    {
        CIquan <- qnorm(1 - (1 - level)/2)
        CIlower <- predVec - CIquan * seVec
        CIupper <- predVec + CIquan * seVec
      
        retMat <- as.matrix(cbind(predVec, CIlower, CIupper))
        colnames(retMat) <- c("Prediction", "Lower CI", "Upper CI")      
        return(retMat)    
    }
    if (interval == "prediction")
    {
        CIquan <- qnorm(1 - (1 - level)/2)
        seVec <- sqrt(object$"gof"[4] + seVec^2)
        CIlower <- predVec - CIquan * seVec
        CIupper <- predVec + CIquan * seVec
      
        retMat <- as.matrix(cbind(predVec, CIlower, CIupper))
        colnames(retMat) <- c("Prediction", "Lower PI", "Upper PI")      
        return(retMat)    
    }
    
    if (pava)
    {
        if (coef(object$object1)[1] > 0)  # works only for a single curve    
        {
            -pava(-predVec)
        } else {
            pava(predVec)
        }
    } else {
        predVec
    }
}  
  

## Calculating a single ED value with confidence interval
"inverseRegBasic" <- function(object, perc, level, interval, method, cgridsize, gridsize, 
lowerRef, upperRef, intType, minmax)
{
#    method <- match.arg(method)
#    intType <- "confidence"  # "prediction"
    predVec <- predict(object)
    if (identical(minmax, "response"))
    {
        maxVal <- max(predVec)   # more sophisticated?
        minVal <- min(predVec)
    } else { # dose-based (better for hormesis data)
        maxVal <- predVec[which.min(object$"dose")]    
        minVal <- predVec[which.max(object$"dose")]    
    }

    ## Using specified lower and upper limits instead of the limits of the response
    if (!is.null(lowerRef)) 
    {
        minVal <- lowerRef
    }
    if (!is.null(upperRef)) 
    {
        maxVal <- upperRef
    }
   
    ## Swapping in case of a decreasing curve (only works for a single curve)
    if (coef(object$"object1")[1] > 0) 
    {
        newPerc <- 100 - perc
    } else {
        newPerc <- perc
    }
    
    ## Scaling percentage up to original response scale
    val <- switch(object$object1$"type",
    "binomial" = newPerc/100,
    "continuous" = (newPerc/100) * (maxVal - minVal) + minVal)    

    ## Checking that ED level can be estimated
    if ( (val < minVal) || (val > maxVal) )
    {
        warning(paste("ED", perc, " cannot be estimated (not sufficient data)", sep = ""), call. = FALSE)
        return(rep(NA, 3))
    }

    ## Searching an initial crude grid
    x <- object$"dose"
    minx <- min(x)
    maxx <- max(x)
    doseDF1 <- data.frame(seq(minx, maxx, length.out = cgridsize))
    colnames(doseDF1) <- as.character(formula(object$object1)[[3]])
    predVec1 <- predict(object, newdata = doseDF1, interval = intType, level = level)
    min1 <- which.min(abs(predVec1[, 1] - val))[1]
    if (identical(interval, "approximate"))
    {
        min1l <- which.min(abs(predVec1[, 2] - val))[1]
        min1u <- which.min(abs(predVec1[, 3] - val))[1]    
    } else {
        min1l <- NA
        min1u <- NA
    }
     
    ## Searching a finer grid
    retVec <- switch(method,
    "bisection" = 
    {  # bisection method
        bisecFct <- function(min1, column)
        {
            if (is.na(min1)) {return(c(NA))}
        
            pFct <- function(x) 
            {
                doseDF <- data.frame(x)
                colnames(doseDF) <- as.character(formula(object$object1)[[3]])
                as.vector(predict(object, newdata = doseDF, interval = intType, 
                level = level)[, column]) - val
            }
            interVal <- c(doseDF1[max(c(min1 - 1, 1)), 1], doseDF1[min(c(min1 + 1, cgridsize)), 1])
#            print(pFct(interVal[2]/2))
            rootVal <- try(uniroot(pFct, interVal)$root, silent = TRUE)
            if (inherits(rootVal, "try-error"))
            {
                warnText <- switch(column,
                "1" = "ED",
                "2" = "Complete confidence interval for ED",  # more informative solution?
                "3" = "Complete confidence interval for ED")
                warning(paste(warnText, perc, " cannot be estimated (bisection method failed)", sep = ""), 
                call. = FALSE)
                return(NA)
            } else {
                return(rootVal)
            }
        } 
        limVec <- c(bisecFct(min1l, 2), bisecFct(min1u, 3))
        if (!any(is.na(limVec))) {limVec <- sort(limVec)}
        c(bisecFct(min1, 1), limVec)        
#        c(bisecFct(min1, 1), sort(c(bisecFct(min1l, 2), bisecFct(min1u, 3)), na.last = TRUE))
    },
    "grid" = 
    {  # grid search (slower, but more robust close to the boundaries of the dose range
        gridSearch <- function(min1, column)
        { 
            if (is.na(min1)) {return(c(NA))}
                    
            doseDF2 <- data.frame(seq(doseDF1[max(c(min1 - 1, 1)), 1], doseDF1[min(c(min1 + 1, cgridsize)), 1], 
            length.out = gridsize))
            colnames(doseDF2) <- as.character(formula(object$object1)[[3]])    
            doseDF2[which.min(abs(predict(object, newdata = doseDF2, interval = intType, 
            level = level)[, column] - val)), 1]
        }
        limVec <- c(gridSearch(min1u, 3), gridSearch(min1l, 2))
        if (!any(is.na(limVec))) {limVec <- sort(limVec)}
        c(gridSearch(min1, 1), limVec)                
#        c(gridSearch(min1, 1), sort(c(gridSearch(min1u, 3), gridSearch(min1l, 2)), na.last = TRUE))
    })
    ## Adjusting confidence limits in case they are on the boundaries of the dose range used
#    if (any(is.na(retVec))) {return(rep(NA, 3))}
    if (is.na(retVec[1])) {return(rep(NA, 3))}  # no limits returned in case no estimate is available
    if ( (!is.na(retVec[2])) && (abs(retVec[2] - minx) < 0.001) ) {retVec[2] <- 0}  # truncation to 0
    if ( (!is.na(retVec[3])) && (abs(retVec[3] - maxx) < 0.001) ) {retVec[3] <- Inf}  # unbounded upper limit    
    
    retVec
}
## Calculating several ED values
"inverseReg" <- Vectorize(inverseRegBasic, "perc")

## Displaying the ED values
EDprint <- function(EDmat, ci, ciLabel, display)
{
    if (display)
    {
        cat("\n")
        cat("Estimated effective doses\n")
        if (!(ci == "none")) 
        {
            ciText <- paste("(", ciLabel, "-based confidence interval(s))\n", sep = "")
            cat(ciText)
        }
        cat("\n") 
        printCoefmat(EDmat)
    }
    invisible(EDmat)  
}

## Calculating ED values for model-robust fits (only works for one curve)
ED.mrdrc <- function(object, respLev, interval = c("none", "approximate", "bootstrap"), level = 0.95, 
method = c("bisection", "grid"), cgridsize = 20, gridsize = 100, display = TRUE, lower = NULL, upper = NULL, 
intType = c("confidence", "prediction"), minmax = c("response", "dose"), n = 1000, seedVal = 200810311, ...)
{
    interval <- match.arg(interval)
    intType <- match.arg(intType)
    method <- match.arg(method)
    minmax <- match.arg(minmax)
  
    ## Calculating by means of inverse regression
    if (identical(interval, "bootstrap"))
    {
        EDvec <- t(inverseReg(object, respLev, level, "none", method, cgridsize, gridsize, lower, upper, intType, 
        minmax))[, 1]
        
        EDci <- EDboot(n, object, respLev, seedVal, level)
        EDmat <- as.matrix(cbind(EDvec, EDci))
    } else {
        EDmat <- t(inverseReg(object, respLev, level, interval, method, cgridsize, gridsize, lower, upper, intType, 
        minmax))
    }
    rownames(EDmat) <- respLev
    
    if (identical(interval, "none"))
    {
        EDmat <- EDmat[, 1, drop = FALSE]  # removing NAs produced by inverseReg()
        colnames(EDmat) <- c("Estimate")
    } else {
        colnames(EDmat) <- c("Estimate", "Lower", "Upper")
        ciLabel <- ifelse(identical(interval, "approximate"), "Approximate variance formula", "Bootstrap")
    }
    EDprint(EDmat, interval, ciLabel, display)
}


plot.mrdrc <- function(x, ..., pava = FALSE)
{
    object <- x
    drcObject <- x$"object1"

    plotPoints.mr <- function(doseVec)
    {
        doseDF <- data.frame(doseVec)
        colnames(doseDF) <- as.character(formula(object$object1)[[3]])  # step needed for loess prediction
        cbind(predict(object, newdata = doseDF, pava = pava))
    }
    drcObject$"curve"[[1]] <- plotPoints.mr
    drcObject$"curve"$"naPlot" <- TRUE  # specifying that no extrapolation will go on when plotting

    plot(drcObject, ...)
}

print.mrdrc <- function(x, ...)
{
    object <- x
    gofVec <- object$"gof"

    cat(paste("\n", "A model-robust dose-response fit", "\n", sep = ""))
    cat(paste("mixing a fit of class '", class(object$"object1")[1], "'",
    " and a fit of class '", class(object$"object2")[1], "'\n\n", sep = ""))
    
    optimalLambda <- object$"lambda"
    cat("Mixing coefficient:", optimalLambda, "\n")
    if (optimalLambda < 1e-5) {noMixing <- TRUE} else {noMixing <- FALSE}
    cat("Degrees of freedom:", object$"df", "\n")
    if (!noMixing)
    {
        cat(paste("(for purely parametric fit: ", format(df.residual(object$"object1"), digits = 5), ")\n\n", sep = ""))
    }
    
    if (identical(object$object1$"type", "binomial"))
    {
        cat("Pearson's chi-square:", format(gofVec[1], digits = 5), "\n")
        if (!noMixing)
        {                
            cat(paste("(for purely parametric fit: ", format(gofVec[2], digits = 5), " (p=", 
            format(1-pchisq(gofVec[2], df.residual(object$"object1")), digits = 2) , "))\n\n", sep = ""))   
        }
        sigma2 <- 0
    } else {
        cat("Residual sum of squares:", format(gofVec[1], digits = 5), "\n")    
        if (!noMixing)
        {
            cat(paste("(for purely parametric fit: ", format(gofVec[2], digits = 5), ")\n\n", sep = ""))
        }
        cat("Residual standard error:", format(gofVec[4], digits = 5), "\n\n")
        
        sigma2 <- 2  # adjustment for AIC ("removing" sigma^2 by subtracting 2)
    } 
    
    cat("AIC:", gofVec[3], "\n")
    if (!noMixing)
    {
        cat(paste("(for purely parametric fit: ", format(AIC(object$"object1") - sigma2, digits = 5), ")\n\n", sep = ""))    
    } else {
        cat("\n")
    }
    invisible(object)
}

"EDboot" <- function(n, object, respLev, seedVal, level)
{
    set.seed(seedVal)

    doseVec <- object$"dose"
    drcObj <- object$object1
    fitVal <- fitted(object)
    lenData <- length(fitVal)
    respType <- drcObj$"type"
    weightsVec <- drcObj$"weights"

    lenED <- length(respLev)

    ## Generating datasets
    dataMat <- switch(respType, 
    "binomial" = t(t(rbinom(n * lenData, weightsVec, fitVal)) / weightsVec),
    "continuous" = rnorm(n * lenData, fitVal, sqrt(object$"gof"[4])))
    dataMat <- matrix(dataMat, nrow = n, byrow = TRUE)
     
    ## Defining apply() functions
    rowFct1 <- function(yVec)
    {
        m1 <- try(drm(yVec ~ doseVec, weights = weightsVec, fct = drcObj$fct, type = respType, 
        start = coef(drcObj)), silent = TRUE)
         
        if (inherits(m1, "try-error"))
        {
            return(rep(NA, lenED))
        }      
        m2 <- loess(yVec ~ doseVec, degree = 1)
        mr <- mrdrm(m1, m2)
         
        as.vector(ED(mr, respLev, display = FALSE))
    }
    rowFct2 <- function(edVec)
    {
        as.vector(quantile(edVec, c((1 - level)/2, 1 - ((1 - level)/2)), na.rm = TRUE))     
    }
    ## Calculating ED values
    if (lenED < 2)
    {
        EDmat <- matrix(rowFct2(apply(dataMat, 1, rowFct1)), nrow = 1)
    } else {
        EDmat <- t(apply(apply(dataMat, 1, rowFct1), 1, rowFct2))
    }        
    EDmat     
}


#mv <- function(n, edval = 50, seedVal = 200810151)
#{
#    set.seed(seedVal)
#
#    doseVec <- deguelin$dose
#    nVec <- deguelin$n
#    if (FALSE)
#    {
##        sigmaVal <- exp.x.mr$
##        lenData <- nrow()
##        doseVec <- exp.x$conc       
##        fittedVec <- fitted(exp.x.mr)          
#    }
#    edVec <- rep(NA, n)
#    lVec <- rep(NA, n)
#    for (i in 1:n)
#    {
#        if (TRUE)
#        {
#            rVec <- rbinom(length(nVec), nVec, fitted(deguelin.mr))
#            m1 <- try(drm(rVec/nVec ~ doseVec, weights = nVec, fct = LL.2(), type = "binomial", start = coef(deguelin.m1)), 
#            silent = TRUE)
#        }
#        if (FALSE)
#        {
#            rVec <- rnorm(lenData, fittedVec, sigmaVal)
#            m1 <- try(drm(rVec ~ doseVec, fct = LL.3(), start = coef(exp.x.m1)), 
#            silent = TRUE)
#        
#        }
#        if (inherits(m1, "try-error"))
#        {
#            edVec[i] <- NA
#        }  else {
#            m2 <- loess(rVec/nVec ~ doseVec, degree = 1)
#            mr <- mrdrm(m1, m2)
#            edVec[i] <- ED(mr, edval, display = FALSE)[1]  
#            lVec[i] <- mr$lambda
#        }
#    }
#    list(edVal = edVec, lambdaVal = lVec)
#}


##  From the R-help e-mail by Ted Harding: http://tolstoy.newcastle.edu.au/R/e2/help/07/03/12853.html
##  See also http://tolstoy.newcastle.edu.au/R/help/05/05/4254.html
"pava" <- function(x, wt = rep(1, length(x)))
{
    n <- length(x)
    if (n <= 1) return(x)
    lvlsets <- 1:n
    repeat 
    {
        viol <- (as.vector(diff(x)) < 0)
        if (!(any(viol))) break
        i <- min((1:(n-1))[viol])

        lvl1 <- lvlsets[i]
        lvl2 <- lvlsets[i+1]
        ilvl <- ( (lvlsets == lvl1) | (lvlsets == lvl2) )

        x[ilvl] <- sum(x[ilvl] * wt[ilvl]) / sum(wt[ilvl])     
        lvlsets[ilvl] <- lvl1
    }
    x
} 


## Examples
if (FALSE)
{
    deguelin.mr <- mrdrm(deguelin.m1, deguelin.m2)
    predict(deguelin.mr, interval="confidence")
    plot(deguelin.m1, ylim=c(0,1))
    lines(deguelin$dose, predict(deguelin.mr), lty=2)

    ## With fixed lambda
    deguelin.mr2 <- mrdrm(deguelin.m1, deguelin.m2, lambda=0.95)
    lines(deguelin$dose, predict(deguelin.mr2), lty=3)

    exp.a.mr <- mrdrm(exp.a.m1, exp.a.m2)
    predict(exp.a.mr, se.fit = TRUE)
    ED(exp.a.mr, c(10, 50, 90))

    ryegrass.mr <- mrdrm(ryegrass.m1, ryegrass.loess)
    predict(ryegrass.mr)
}


#press(deguelin.m1, deguelin.m2, seq(0, 1, length.out=10), criterion="ll")


## trace hat component in drm() model fit? 
#traceHat.drc <- function(object)
#{
#    object$"sumList"$lenData - df.residual(object)
#}
#
#traceHat.loess <- function(object)
#{
#     object$trace.hat
#}




## Example R lines
if (FALSE)
{

deguelin.m1<-drm(r/n~dose, weights=n, data=deguelin, fct=LL.2(), type="binomial")

deguelin.m2<-loess(r/n~dose, data=deguelin, degree=1)
vPress(deguelin.m1, deguelin.m2, seq(0, 1, length.out=10))

deguelin.m3<-loess(r/n~log10(dose), data=deguelin, degree=1)
vPress(deguelin.m1, deguelin.m3, seq(0, 1, length.out=10))

}

## Using lmeSplines
if (FALSE)
{
library(lmeSplines)

## lettuce
let2<-lettuce
let2<-cbind(let2, all=rep(1, 14))
let2$Zt<-smspline(~conc, data=lettuce)

let2.m1<-lme(weight~conc, data=let2, random=list(all=pdIdent(~Zt-1)))
plot(weight~conc, data=lettuce, log="x")
lines(let2$conc, fitted(let2.m1))


## deguelin
deg2<-deguelin
deg2<-cbind(deg2, all=rep(1, 6))
deg2$Zt<-smspline(~dose, data=deg2)

deg2.m1<-lme(r/n~dose, data=deg2, random=list(all=pdIdent(~Zt-1)))
plot(r/n~dose, data=deg2)
lines(deg2$dose, fitted(deg2.m1))

## bin.mat
bm2<-bin.mat[c(3,6,9,12,15),]
bm2<-cbind(bm2, all=rep(1, 5))
bm2$Zt<-smspline(~conc, data=bm2)

bm2.m1<-lme(matured/total~conc, data=bm2, random=list(all=pdIdent(~Zt-1)))
plot(matured/total~conc, data=bm2)
lines(bm2$conc, fitted(bm2.m1))  # straight line!



## exp.a
ea2<-exp.a
ea2<-cbind(ea2, all=rep(1,54))

ea2$Zt<-smspline(~x, data=ea2)
ea2.m1<-lme(y~x, data=ea2, random=list(all=pdIdent(~Zt-1)))
plot(y~x, ea2, log="x")
lines(ea2$x, fitted(ea2.m1))
summary(ea2.m1)


ea2.m2<-drm(y~x, data=ea2, fct=LL.3())
plot(ea2.m2, add=T, lty=2)

cVec<-ea2$x[-c(7:12)]
Zt2<-smspline(cVec)
ea2$Zt2<-approx.Z(Zt2, cVec, ea2$x)
ea2.m3<-lme(y~x, data=ea2, random=list(all=pdIdent(~Zt2-1)))
lines(ea2$x, fitted(ea2.m3), col=2)

cVec<-ea2$x[-c(13:18)]
Zt2<-smspline(cVec)
ea2$Zt2<-approx.Z(Zt2, cVec, ea2$x)
ea2.m3<-lme(y~x, data=ea2, random=list(all=pdIdent(~Zt2-1)))
lines(ea2$x, fitted(ea2.m3), col=5)



}