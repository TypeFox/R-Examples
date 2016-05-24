### ===== actuar: An R Package for Actuarial Science =====
###
### Sample empirical limited value functions for individual and
### grouped data.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca> and
###          Mathieu Pigeon

elev <- function(x, ...)
{
    Call <- match.call()
    UseMethod("elev")
}

elev.default <- function(x, ...)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()
    chkDots(...)                        # method does not use '...'
    FUN <- function(limit)
        sapply(limit, function(x, y) mean(pmin(x, y)), x = x)
    environment(FUN) <- new.env()
    assign("x", sort(x), envir = environment(FUN))
    assign("n", length(unique(x)), envir = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    attr(FUN, "call") <- Call
    attr(FUN, "grouped") <- FALSE
    FUN
}

### Function 'elev.grouped.data' below returns a function that uses
### data stored in its environment. Avoid false positive in R CMD
### check.
if(getRversion() >= "2.15.1")  utils::globalVariables(c("cj", "nj"))

### This function assumes right-closed intervals, but the numerical
### values are identical for left-closed intervals.
elev.grouped.data <- function(x, ...)
{
    if (!exists("Call", inherits = FALSE))
        Call <- match.call()
    chkDots(...)                        # method does not use '...'
    FUN <- function(limit)
    {
        ## Explicitely get the data from the function environment.
        ## cj <- eval(expression(cj))
        ## nj <- eval(expression(nj))

        ## Number of classes.
        r <- length(nj)

        ## This is to avoid numerical problems.
        limit <-  pmin(limit, cj[r + 1L])

        ## Class in which the limit is located.
        cl <- findInterval(limit, cj, all.inside = TRUE)

        ## Means for all classes below each limit.
        cjt <- head(cj, max(cl))        # upper bounds
        res1 <- sapply(cl - 1L, function(n, x)
                       drop(crossprod(head(x, n), head(nj, n))),
                       (head(cjt, -1) + tail(cjt, -1))/2)

        ## Means for classes with each limit.
        cjt <- cj[cl]                   # lower bounds
        njt <- nj[cl]                   # frequencies
        p <- (limit - cjt) / (cj[cl + 1L] - cjt) # prop. to take
        res2 <- njt * p * (cjt + limit)/2 + njt * (1 - p) * limit

        ## Means for classes above each limit.
        res3 <- limit * sapply(r - cl, function(n, x) sum(tail(x, n)),
                               tail(nj, -min(cl)))

        ## Total
        (res1 + res2 + res3)/sum(nj)
    }

    environment(FUN) <- new.env()
    assign("cj", eval(expression(cj), envir = environment(x)),
           envir = environment(FUN))
    assign("nj", x[, 2L], envir = environment(FUN))
    assign("n", nrow(x), envir = environment(FUN))
    class(FUN) <- c("elev", class(FUN))
    attr(FUN, "call") <- Call
    attr(FUN, "grouped") <- TRUE
    FUN
}

### Essentially identical to stats::print.ecdf().
print.elev <- function(x, digits = getOption("digits") - 2, ...)
{
    ## Utility function
    numform <- function(x) paste(formatC(x, digits = digits), collapse = ", ")

    ## The rest is adapted from ecdf()
    varname <- if (attr(x, "grouped")) "cj" else "x"
    cat("Empirical LEV \nCall: ")
    print(attr(x, "call"), ...)
    n <- length(xx <- eval(parse(text = varname), envir = environment(x)))
    i1 <- 1L:min(3L, n)
    i2 <- if (n >= 4L) max(4L, n - 1L):n else integer(0)
    cat(" ", varname, "[1:", n, "] = ", numform(xx[i1]), if (n > 3L) ", ",
        if (n > 5L) " ..., ",  numform(xx[i2]), "\n", sep = "")
    invisible(x)
}

### Essentially identical to stats::summary.ecdf().
summary.elev <- function (object, ...)
{
    cat("Empirical LEV:\t ", eval(expression(n), envir = environment(object)),
        "unique values with summary\n")
    summary(knots(object), ...)
}

### Essentially identical to stats::knots.stepfun().
knots.elev <- function(Fn, ...)
{
    if (attr(Fn, "grouped"))
        eval(expression(cj), envir = environment(Fn))
    else
        eval(expression(x), envir = environment(Fn))
}

plot.elev <- function(x, ..., main = NULL, xlab = "x", ylab = "Empirical LEV")
{
    if (missing(main))
        main <- {
            cl <- attr(x, "call")
            deparse(if (!is.null(cl)) cl else sys.call())
        }

    kn <- knots(x)
    Fn <- x(kn)
    plot(kn, Fn,  ..., main = main, xlab = xlab, ylab = ylab)
}
