"confint.drc" <- function(object, parm, level = 0.95, pool = TRUE, ...)
#"confint.drc" <- function(object, parm, level = 0.95, type = "t", pool = TRUE, ...)
{
    ## Matching parameter names
    if (!missing(parm))
    {
        matchVec <- object$"parNames"[[2]] %in% parm
        if (!any(matchVec)) {stop("The 'parm' argument does not match an actual parameter")}
    } else {
        matchVec <- rep(TRUE, length(object$"parNames"[[2]]))
    }

    ## Constructing matrix of confidence intervals    
    confint.basic(summary(object, pool = pool)$"coefficients"[matchVec, 1:2, drop = FALSE], 
    level, object$"type", df.residual(object))  
    
#    ## Retrieving estimates and estimated standard errors
#    estMat <- summary(object, pool = pool)$"coefficients"[matchVec, 1:2, drop = FALSE]
#
#    ## Constructing matrix of confidence intervals
#    confMat <- matrix(0, dim(estMat)[1], 2)
#
#    alphah <- (1 - level)/2 
#    if (type == "u") {two <- qnorm(1 - alphah)}
#    if (type == "t") {two <- qt(1 - alphah, df.residual(object))}
#    confMat[, 1] <- estMat[, 1] -  two * estMat[, 2]
#    confMat[, 2] <- estMat[, 1] +  two * estMat[, 2]
#
#    ## Formatting matrix
#    colnames(confMat) <- c(paste(format(100 * alphah), "%", sep = " "), paste(format(100*(1 - alphah)), "%", sep = " ") )
#    rownames(confMat) <- rownames(estMat)
#
#    return(confMat)  
}

## Defining basic function for providing confidence intervals
"confint.basic" <- function(estMat, level, intType, dfres, formatting = TRUE)
{
    alphah <- (1 - level)/2 
#    if (type == "u") {two <- qnorm(1 - alphah)}
#    if (type == "t") {two <- qt(1 - alphah, df.residual(object))}    
    tailPercentile <- switch(intType, 
    binomial = qnorm(1 - alphah), 
    continuous = qt(1 - alphah, dfres),
    event = qnorm(1 - alphah),
    Poisson = qnorm(1 - alphah))
    
    estVec <- estMat[, 1]
    halfLength <- tailPercentile * estMat[, 2]
    confMat <- matrix(c(estVec -  halfLength, estVec +  halfLength), ncol = 2)    
    
    ## Formatting matrix
    if (formatting)
    {
        colnames(confMat) <- c(paste(format(100 * alphah), "%", sep = " "), paste(format(100*(1 - alphah)), "%", sep = " "))
        rownames(confMat) <- rownames(estMat)
    }

    return(confMat)    
}
