"diagnostics" <- function(object)
{
    ## Displaying information about the estimation
    if (!inherits(object, "drc")) {stop("Object not of class 'drc'")}

    convStatus <- object$"fit"$convergence
    if (convStatus)  # (convStatus == 0)
    {
        convTxt <- "Estimation procedure converged"
    } else {
        convTxt <- "Estimation procedure failed to converged"
    }
    cat(paste(convTxt, ".\n\n", sep=""))

    convRate <- object$"fit"$counts
    cat("Convergence measured in number of evaluations:\n\n")

    cat(paste("Function  :", convRate[1], "\n"))
    cat(paste("Derivative:", convRate[2], "\n\n"))

    options(warn=-1)
    if (any(is.na(summary(object)$"estimates"[,2]))) {cat("Some estimated standard errors not available.\n\n")}
    options(warn=0)
}
