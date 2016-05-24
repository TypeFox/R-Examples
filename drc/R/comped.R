"comped" <- function(est, se, log = TRUE, interval = TRUE, operator = c("-", "/"), level = 0.95, df = NULL)
{
    operator <- match.arg(operator)
    if (identical(operator, "-")) 
    {
        opText <- "difference"
    } else {
        opText <- "ratio"
    }
    
    vcMat <- diag(se^2)
    if (!log)
    {
        if (identical(operator, "-"))
        {
#            resList <- compute.delta.method(vcMat, expression(b1-b2), est, c("b1", "b2"), print = FALSE)
            
            derivVec <- c(1, -1)
            estVal <- est[1] - est[2]
        } else { 
#            resList <- compute.delta.method(vcMat, expression(b1/b2), est, c("b1", "b2"), print = FALSE)
            
            derivVec <- c(1 / est[2], -est[1] / (est[2]^2))           
            estVal <- est[1] / est[2]
        }
    } else {
 #       resList <- compute.delta.method(vcMat, expression(b1-b2), est, c("b1", "b2"), print = FALSE)
 
        derivVec <- c(1, -1)
        estVal <- est[1] - est[2]
    }
    resList <- list(estimate = estVal, se = sqrt(as.numeric(derivVec %*% vcMat %*% derivVec)))            
    colNames <- c("Estimate", "Std. Error")    
    
    edMat <- matrix(c(resList$estimate, resList$se), nrow = 1)
    if (interval)
    {
        colNames <- c(colNames, "Lower", "Upper")
        
        ## Setting degrees of freedom (by default based on normality)
        if (is.null(df)) 
        {
            quanVal <- qnorm(1 - (1 - level)/2)
        } else {
            quanVal <- qt(1 - (1 - level)/2, df)
        }
        ciMat <- matrix(c(edMat[1, 1] - quanVal * edMat[1, 2], edMat[1, 1] + quanVal * edMat[1, 2]), nrow = 1)        
        
        if (log && (identical(operator, "/")))
        {
            ciMat <- exp(ciMat)
        }
        edMat <- cbind(edMat, ciMat)
    }
    colnames(edMat) <- colNames

    cat("\n")
    cat("Estimated", opText, "of effective doses\n")
    if (interval && log && (identical(operator, "/"))) {cat("(confidence interval on original scale)\n")}
    cat("\n") 

    printCoefmat(edMat)
    invisible(edMat)
}