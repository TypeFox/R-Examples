"print.summary.drc" <-
function(x, ...)
{
    object <- x
    
    cat("\n")
#    cat(paste("Model fitted: ", object$"fctName", "\n", sep = ""))

    if (!is.null(object$"noParm"))
    {
        modelText <- paste("Model fitted: ", object$"text", " (", object$"noParm", " parms)", "\n", sep = "")
    } else {
        modelText <- paste("Model fitted: ", object$"text", "\n", sep = "")
    }
    cat(modelText)

    if (!is.null(object$"robust"))
    {
        cat("\n")
        cat("Robust estimation:", object$"robust", "\n")   
    }

    cat("\n")
    cat("Parameter estimates:\n\n")
    printCoefmat(object$"coefficients")

    if (!is.na(object$"resVar"))
#    if ((!is.null(object$"resVar")) && (!identical(object$"type", "binomial")))
    {
        cat("\nResidual standard error")
        
        rseMat <- object$"rseMat"
        if (nrow(rseMat) > 1)
        {
            cat("s:\n\n")
        } else {
            cat(":\n\n")
        }

        printFct <- function(x)
        {
            paste(format(x[1]), paste("(", as.character(x[2]), " degrees of freedom)\n", sep = ""))
        }
        cat(paste(rownames(rseMat), apply(rseMat, 1, printFct)), sep = "")
    }

    if (!is.null(object$"varComp"))
    {
        cat("\n")
        cat("Estimated variance components:\n\n")
        printCoefmat(object$"varComp")
    }
    
    if (!is.null(object$"varParm"))  # summary of variance as power of mean
    {
        if (object$"varParm"$"type" == "varPower")
        {
            cat("\n")
            cat("Heterogeneity adjustment through power-of-the-mean variance model\n\n")    
            if (dim(object$"varParm"$"estimates")[1] > 1)
            {
                printCoefmat(object$"varParm"$"estimates"[2, , drop=FALSE])  
                # only displaying power exponent, not residual variance
            } else {
                printCoefmat(object$"varParm"$"estimates")
            }
        }
        if (object$"varParm"$"type" == "hetvar")
        {
            cat("\n")
            cat("Estimated heterogeneous variances:\n\n")    
            printCoefmat(object$"varParm"$"estimates"[, 1, drop = FALSE])
        }        
    }
    
    if (!is.null(object$"boxcox"))  # summary of Box-Cox transformation
    {
        bcObj <- object$"boxcox"
        lambda <- bcObj$"lambda"
        if (is.na(lambda)) 
        {
            # empty
        } else {
#            pVal <- format(object$"boxcox"[2], digits=3)
#            boxcoxci <- c(format(ci[1], digits = 3), format(ci[2], digits = 3))

            cat("\n")
            cat("Non-normality/heterogeneity adjustment through optimal Box-Cox transformation\n\n")

            ci <- bcObj$"ci"
            if (!is.na(ci[1]))
            {        
                cat("Estimated lambda:", format(lambda, digits = 3), "\n")
#                cat("P-value for test of null hypothesis that lambda=1:", pVal, "\n")
                ci <- format(ci, digits = 3) 
                ciStr <- paste("[", ci[1], ",", ci[2], "]", sep="")
                cat("Confidence interval for lambda:", ciStr, "\n\n")
            } else {
                cat("Specified lambda:", lambda, "\n\n")        
            }
        }
    }
         
    invisible(object)
}
