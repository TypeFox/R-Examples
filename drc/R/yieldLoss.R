"yieldLoss" <- function(object, interval = c("none", "as"), 
level = 0.95, display = TRUE)
{
    interval <- match.arg(interval)

    fixedFct <- genFixedFct(object$"fct"$"fixed")
    indMat <- object$"indexMat"
    parmMat <- object$"parmMat"
    rowNames <- colnames(parmMat)
    vcMat <- vcov(object)
    
    ## deriv(~d, c("b", "c", "d", "e", "f"), function(b,c,d,e,f){}   
    Afct <- function(parm)
    {
        .value <- fixedFct(parm)[3]
        .grad <- array(0, c(length(.value), 5L), list(NULL, c("b", "c", "d", "e", "f")))
        .grad[, "b"] <- 0
        .grad[, "c"] <- 0
        .grad[, "d"] <- 1
        .grad[, "e"] <- 0
        .grad[, "f"] <- 0
        attr(.value, "gradient") <- .grad
        .value
    }
    
    ## deriv(~d/e, c("b", "c", "d", "e", "f"), function(b,c,d,e,f){})
    Ifct <- function(parm) 
    {
        parm2 <- fixedFct(parm)    
        d <- parm2[3]
        e <- parm2[4]

        .value <- d/e
        .grad <- array(0, c(length(.value), 5L), list(NULL, c("b", "c", "d", "e", "f")))
        .grad[, "b"] <- 0
        .grad[, "c"] <- 0
        .grad[, "d"] <- 1/e
        .grad[, "e"] <- -(d/e^2)
        .grad[, "f"] <- 0
        attr(.value, "gradient") <- .grad
        .value
    }

    ncim <- ncol(indMat)
    Amat <- matrix(NA, ncim, 2) 
    Imat <- matrix(NA, ncim, 2)     
    for (i in 1:ncim)
    {
        indVec <- indMat[, i]
        parmVeci <- parmMat[, i]
        vcMati <- vcMat[indVec, indVec]
        
        aDerVal <- Afct(parmVeci)
        Amat[i, 1] <- aDerVal
        aDer <- fixedFct(attr(aDerVal, "gradient"), allComp = FALSE) 
        Amat[i, 2] <- sqrt(aDer %*% vcMati %*% aDer)
        
        IDerVal <- Ifct(parmVeci)
        Imat[i, 1] <- IDerVal
        IDer <- fixedFct(attr(IDerVal, "gradient"), allComp = FALSE) 
        Imat[i, 2] <- sqrt(IDer %*% vcMati %*% IDer)        
    }
    
    if (identical(interval, "none"))
    {
        colNames <- c("Estimate", "Std. Error")
        intervalLabel <- NULL       
    }
    if (identical(interval, "as"))
    {        
        ciFct <- function(resMat)
        {
            quanVal <- ifelse(identical(object$"type", "continuous"), 
            qt(1 - (1 - level)/2, df.residual(object)), qnorm(1 - (1 - level)/2))
            
            rm1 <- resMat[, 1]
            rm2 <- resMat[, 2]                            
            ciMat <- matrix(0, ncim, 2)
            ciMat[, 1] <- rm1 - quanVal * rm2
            ciMat[, 2] <- rm1 + quanVal * rm2

            as.matrix(cbind(resMat, ciMat))
        }
        Amat <- ciFct(Amat)
        Imat <- ciFct(Imat)
        
        colNames <- c("Estimate", "Std. Error", "Lower", "Upper")
        intervalLabel <- "Asymptotically"
    }    
    
    ## Assigning column and row names
    colnames(Amat) <- colNames
    colnames(Imat) <- colNames    
    rownames(Amat) <- rowNames
    rownames(Imat) <- rowNames
        
    ## Displaying and returning output    
    resPrint(Amat, "Estimated A parameters", interval, intervalLabel, display = display)
    resPrint(Imat, "\nEstimated I parameters", interval, intervalLabel, display = display)        
    invisible(list(A = Amat, I = Imat))
}

genFixedFct <- function(fixed)
{
    notFixed <- is.na(fixed)
    numParm <- length(fixed)
    
    function(parm, allComp = TRUE)
    {
        if (allComp)
        {
            parmVec <- rep(NA, numParm)
            parmVec[!notFixed] <- fixed[!notFixed]
            parmVec[notFixed] <- parm
            parmVec
        } else {
            parm[notFixed] 
        }
    }
}