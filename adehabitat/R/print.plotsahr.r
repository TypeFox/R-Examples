"print.plotsahr" <- function(x, ...)
{
    cat("***** List of class \"plotsahr\" *****\n\n")
    cat("Selection ratios are computed for the following variables:\n\n")
    for (i in 1:length(x))
        cat(names(x)[i], "\n")
    cat("each variable is a component of the list\n\n")
}

