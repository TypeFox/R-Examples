##
## Create the varying part of a fitted value profile.  Used in 'predVals', not
## meant to be called directly by users.
##
## ARGUMENTS:
##   data: data frame containing raw data
##   cols: columns of varying variables
##   xlims: named list specifying limits to use for continuous variables
##     (default is to use observed range)
##   n: number of grid points for continuous variables
##
## RETURN:
##   Data frame containing grid values for specified variables
##
getXvals <- function(data, cols, xlims, n)
{
    ans <- list()
    for (i in seq_along(cols)) {
        xcol <- cols[i]
        xname <- names(data)[xcol]
        x <- data[, xcol]

        ## Construct the sequence of values of the variable to evaluate at
        if (all(unique(x) %in% c(0, 1))) {
            ## Binary variable: use 0 and 1
            xs <- c(0, 1)
        } else if (is.numeric(x)) {
            ## Continuous variable: make a grid of size 'n' within the specified
            ## limits (or observed range if none specified)
            xlim <- xlims[[xname]]
            if (is.null(xlim))
                xlim <- range(x)
            xs <- seq(xlim[1], xlim[2], length.out = n)
        } else if (is.factor(x)) {
            ## Factor variable: use all levels
            xs <- factor(levels(x), levels = levels(x))
        } else if (is.logical(x)) {
            ## Logical variable: use TRUE and FALSE
            xs <- c(FALSE, TRUE)
        }

        ans[[i]] <- xs
        names(ans)[i] <- xname
    }

    ## Create all possible combinations of the specified variables
    ans <- expand.grid(ans)
    return(ans)
}
