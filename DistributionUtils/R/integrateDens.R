integrateDens <- function(densFn = "norm", lower = -Inf, upper = Inf,
                          subdivisions = 100, ...){
    ## Uses code from distCheck by Diethelm Wuertz, modified by David Scott

    CALL <- match.call()
    dfun <- match.fun(paste("d", densFn, sep = ""))

    ## Integrate density over range
    totalPr <- integrate(dfun, lower = lower, upper = upper,
                         subdivisions = subdivisions,
                         stop.on.error = FALSE, ...)

    # Return Value:
    return(totalPr)
}


