"summary.drc" <-
function(object, od = FALSE, pool = TRUE, ...)
{
    ## Producing a summary of a model fit
#    sumVec1 <- object$fit  # object[[2]]
#    sumVec2 <- object$summary  # object[[4]]
#    parNames <- object$"parNames"[[1]]  # object[[6]]

    ## Calculating variance-covariance matrix from Hessian
#    em <- object$"estMethod"
#    parVec <- (em$parmfct)(object$fit, fixed = FALSE)
    parVec <- as.vector(coef(object))
#    notNA <- !is.na(parVec) 
#    varMat <- (object$"scaleFct")( (em$vcovfct)(object) )
    varMat <- vcov(object, od = od, pool = pool)
        
    ## Calculating estimated residual variance 
    ## and unscaled variance-covariance matrix
#    if (!is.null(em$rvfct)) 
#    {
#        resVar <- (em$rvfct)(object)
#    } else {
#        resVar <- NULL
#    }    
    resVar <- rse(object, TRUE)
    if (!is.null(resVar))
    {
        varMat.us <- varMat / (2*resVar)    
    } else {
        varMat.us <- NULL
    }
   
    ## Calculating the residual standard error(s)
    if ((!is.null(object$"objList")) && (!pool))
    {
        objList <- object$"objList"
        lenol <- length(objList)
        
        rseMat <- matrix(NA, lenol, 2)
        rownames(rseMat) <- names(objList) 
        resVar <- as.vector(unlist(lapply(objList, rse, resvar = TRUE)))
        rseMat[, 1] <- sqrt(resVar)  # only to keep resVar
        rseMat[, 2] <- as.vector(unlist(lapply(objList, df.residual)))
    } else {
        resVar <- rse(object, TRUE)
        
        rseMat <- matrix(NA, 1, 2)
        rownames(rseMat) <- ""
        rseMat[1, 1] <- sqrt(resVar)
        rseMat[1, 2] <- df.residual(object)
    }
    colnames(rseMat) <- c("rse", "df")
    
#    ## Adjusting for over-dispersion using the Pearson statistic
#    if (od && (!is.null(object$"gofTest"))) 
#    { 
#        varMat <- varMat*(object$"gofTest"[1]/object$"gofTest"[2])
#    }
    estSE <- sqrt(diag(varMat))

    ## Calculating estimated standard errors for robust methods
    
    ## M-estimators
    if (!is.null(object$robust) && object$robust%in%c("metric trimming", "metric Winsorizing", "Tukey's biweight"))
    {
        psi.trimmed <- function(u, deriv = 0)
        {
            if (deriv == 0)
            {
                retVec <- u
                retVec[ abs(u) > 1.345 ] <- 0
            }
            if (deriv == 1)
            {
                retVec <- rep(1, length(u))
                retVec[ abs(u) > 1.345 ] <- 0
            }
            return(retVec)            
        }
    
    
        if (object$robust=="Tukey's biweight")
        {
            psifct <- psi.bisquare  # in MASS
        }
        if (object$robust=="metric Winsorizing")
        {
            psifct <- psi.huber  # in MASS
        }
        if (object$robust=="metric trimming")
        {
            psifct <- psi.trimmed
        }
    
#        resVec <- residuals(object)
        resVec <- (object)[["predres"]][, "Residuals"]
        psiprime <- psifct(resVec/sqrt(resVar), deriv=1)
        meanpp <- mean(psiprime)
        
        notNA <- !is.na(parVec) 
        sumVec1 <- object$fit 
        K <- 1 + length(parVec[notNA])*var(psiprime)/(object$summary[7]*meanpp^2)        
        w <- psifct(resVec/sqrt(resVar))
        s <- sum((resVec*w)^2)/object$summary[6]
        print(parVec[notNA])
        print(var(psiprime))
        print(c(K,w,s))
        stddev <- sqrt(s)*(K/meanpp)       
        estSE <- sqrt(diag(solve(sumVec1$hessian[notNA, notNA]/mean(psiprime)*resVar)))*stddev
    }


    ## Forming a matrix of results        
#    resultMat <- matrix(0, sum(notNA), 4, dimnames = list(parNames, c("Estimate", "Std. Error", "t-value", "p-value")))    
#    resultMat[, 1] <- parVec[notNA]
    
    parNames <- object$"parNames"[[1]]    
    resultMat <- matrix(NA, length(parVec), 4, 
    dimnames = list(parNames, c("Estimate", "Std. Error", "t-value", "p-value")))    

    resultMat[, 1] <- parVec
    resultMat[, 2] <- estSE
    tempStat <- resultMat[, 1] / resultMat[, 2]
    resultMat[, 3] <- tempStat
    
    ## Using t-distribution for continuous data
    ##  only under the normality assumption
    if (object$"type" == "continuous")
    {
        pFct <- function(x) {pt(x, df.residual(object))}
    } else {
        pFct <- pnorm
    }    
    resultMat[, 4] <- pFct(-abs(tempStat)) + (1 - pFct(abs(tempStat)))

    ## Separating out variance parameters
if (FALSE)
{    
    if (!is.null(object$"varParm"))
    {
        indexVec <- object$"varParm"$"index"
        varParm <- object$"varParm"

        estVec <- resultMat[-indexVec, , drop = FALSE]
        if (object$"varParm"$"type" == "varPower")
        {
            estVec[2, 3] <- (estVec[2, 1] - 0)/estVec[2, 2]  # testing the hypothesis theta=0
            estVec[2, 4] <- 2*pt(-abs(estVec[2, 3]), df.residual(object))
        }
        varParm$"estimates" <- estVec
        
        resultMat <- resultMat[indexVec,]
        varMat <- varMat[indexVec, indexVec]  # for use in ED/MAX/SI
    } else {
        varParm <- NULL
    }
}

    fctName <- deparse(object$call$fct)    

    sumObj <- list(resVar, varMat, resultMat, object$"boxcox", fctName, object$"robust", NULL, object$"type", 
    df.residual(object), varMat.us, object$"fct"$"text", object$"fct"$"noParm", rseMat)
    names(sumObj) <- c("resVar", "varMat", "coefficients", "boxcox", "fctName", "robust", "varParm", "type", 
    "df.residual", "cov.unscaled", "text", "noParm", "rseMat")
    
    class(sumObj) <- c("summary.drc")
    return(sumObj)
}
