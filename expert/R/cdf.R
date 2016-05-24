### ===== expert =====
###
### cdf() returns a function object to compute the cumulative
### distribution function based on the results of expert()
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Mathieu Pigeon <mathieu.pigeon.3@ulaval.ca>

cdf <- function(x, ...)
{
    ## Compute the cumulative distribution function from an
    ## object of class 'expert'.
    if (!inherits(x, "expert"))
        stop("'x' must be an object of class \"expert\"")

    y <- x$probs
    x <- x$breaks

    ## Create an object of class 'cdf'.
    res <- approxfun(x, c(0, pmin(cumsum(y), 1)), yleft = 0, yright = 1,
                     method = "constant", ties = "ordered")
    class(res) <- c("cdf", "stepfun", class(res))
    attr(res, "call") <- sys.call()
    res
}

### Essentially identical to stats::print.ecdf().
print.cdf <- function(x, digits = getOption("digits") - 2, ...)
{
    ## Utility function
    numform <- function(x) paste(formatC(x, dig = digits), collapse = ", ")

    ## The rest is adapted from ecdf()
    cat("Aggregate Expert CDF\nCall: ")
    print(attr(x, "call"), ...)
    nc <- length(xxc <- get("x", env = environment(x)))
    nn <- length(xxn <- get("y", env = environment(x)))
    i1 <- 1:min(3, nc)
    i2 <- if (nc >= 4) max(4, nc - 1):nc else integer(0)
    i3 <- 1:min(3, nn)
    i4 <- if (nn >= 4) max(4, nn - 1):nn else integer(0)
    cat("    x = ", numform(xxc[i1]), if (nc > 3) ", ",
        if (nc > 5) " ..., ", numform(xxc[i2]), "\n", sep = "")
    cat(" F(x) = ", numform(xxn[i3]), if (nn > 3) ", ",
        if (nn > 5) " ..., ", numform(xxn[i4]), "\n", sep = "")
    invisible(x)
}

### Identical to stats::knots.stepfun().
knots.cdf <- stats:::knots.stepfun

### Essentially identical to stats::print.ecdf().
plot.cdf <- function(x, ..., ylab = "F(x)", verticals = FALSE,
                     col.01line = "gray70")
{
    plot.stepfun(x, ..., ylab = ylab, verticals = verticals)
    abline(h = c(0,1), col = col.01line, lty = 2)
}
