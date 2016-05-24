### ===== expert =====
###
### Print method for object of class "expert"
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
###          Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

print.expert <- function(x, ...)
{
    ## Utility function
    numform <- function(x, w)
        formatC(x, digits = 2, width = w, format = "fg")

    breaks <- x$breaks
    probs <- x$probs
    n <- length(breaks)
    alpha <- x$alpha

    ## Display results using a data frame with formatted class in the
    ## first column.
    w <- max(nchar(formatC(breaks[-1], format = "fg"))) # longest upper bound
    fbreaks <- paste("(", numform(breaks[-n], -1),
                     ", ", numform(breaks[-1], w), "]", sep = "")
    res <- data.frame(Interval = fbreaks, Probability = probs)

    cat("\nAggregated Distribution Using", attr(x, "method"), "Model\n\n")
    print(res, row.names = FALSE, ...)
    if (!is.null(alpha))
        cat("\n Alpha:", format(alpha, ...), "\n")
    invisible(x)
}
