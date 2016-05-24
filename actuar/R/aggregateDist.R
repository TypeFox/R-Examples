### ===== actuar: An R Package for Actuarial Science =====
###
### Use one of five methods to compute the aggregate claim amount
### distribution of a portfolio over a period given a frequency and a
### severity model or the true moments of the distribution.
###
### AUTHORS: Vincent Goulet <vincent.goulet@act.ulaval.ca>,
### Louis-Philippe Pouliot

aggregateDist <-
    function(method = c("recursive", "convolution", "normal", "npower", "simulation"),
             model.freq = NULL, model.sev = NULL, p0 = NULL, x.scale = 1,
             convolve = 0, moments, nb.simul, ...,
             tol = 1e-06, maxit = 500, echo = FALSE)
{
    Call <- match.call()

    ## The method used essentially tells which function should be
    ## called for the calculation of the aggregate claims
    ## distribution.
    method <- match.arg(method)

    if (method == "normal")
    {
        ## An error message is issued if the number of moments listed
        ## is not appropriate for the method. However it is the user's
        ## responsability to list the moments in the correct order
        ## since the vector is not required to be named.
        if (missing(moments) || length(moments) < 2)
            stop("'moments' must supply the mean and variance of the distribution")
        FUN <- normal(moments[1], moments[2])
        comment(FUN) <- "Normal approximation"
    }

    else if (method == "npower")
    {
        if (missing(moments) || length(moments) < 3)
            stop("'moments' must supply the mean, variance and skewness of the distribution")
        FUN <- npower(moments[1], moments[2], moments[3])
        comment(FUN) <- "Normal Power approximation"
    }
    else if (method == "simulation")
    {
        if (missing(nb.simul))
            stop("'nb.simul' must supply the number of simulations")
        if (is.null(names(model.freq)) && is.null(names(model.sev)))
            stop("expressions in 'model.freq' and 'model.sev' must be named")
        FUN <- simS(nb.simul, model.freq = model.freq, model.sev = model.sev)
        comment(FUN) <- "Approximation by simulation"
    }
    else
    {
        ## "recursive" and "convolution" cases. Both require a
        ## discrete distribution of claim amounts, that is a vector of
        ## probabilities in argument 'model.sev'.
        if (!is.numeric(model.sev))
            stop("'model.sev' must be a vector of probabilities")

        ## Recursive method uses a model for the frequency distribution.
        if (method == "recursive")
        {
            if (is.null(model.freq) || !is.character(model.freq))
                stop("frequency distribution must be supplied as a character string")
            dist <- match.arg(tolower(model.freq),
                              c("poisson", "geometric", "negative binomial",
                                "binomial", "logarithmic"))
            FUN <- panjer(fx = model.sev, dist = dist, p0 = p0,
                          x.scale = x.scale, ..., convolve = convolve,
                          tol = tol, maxit = maxit, echo = echo)
            comment(FUN) <- "Recursive method approximation"
        }

        ## Convolution method requires a vector of probabilites in
        ## argument 'model.freq'.
        else if (method == "convolution")
        {
            if (!is.numeric(model.freq))
                stop("'model.freq' must be a vector of probabilities")
            FUN <- exact(fx = model.sev, pn = model.freq, x.scale = x.scale)
            comment(FUN) <- "Exact calculation (convolutions)"
        }
        else
            stop("internal error")
    }

    ## Return cumulative distribution function
    class(FUN) <- c("aggregateDist", class(FUN))
    attr(FUN, "call") <- Call
    FUN
}

print.aggregateDist <- function(x, ...)
{
    cat("\nAggregate Claim Amount Distribution\n")
    cat("  ", label <- comment(x), "\n\n", sep = "")

    cat("Call:\n")
    print(attr(x, "call"), ...)
    cat("\n")

    if (label %in% c("Exact calculation (convolutions)",
                     "Recursive method approximation",
                     "Approximation by simulation"))
    {
        n <- length(get("x", envir = environment(x)))
        cat("Data:  (", n, "obs. )\n")
        numform <- function(x) paste(formatC(x, digits = 4, width = 5), collapse = ", ")
        i1 <- 1L:min(3L, n)
        i2 <- if (n >= 4L)
            max(4L, n - 1L):n
        else integer()
        xx <- eval(expression(x), envir = environment(x))
        cat(" x[1:", n, "] = ", numform(xx[i1]), if (n > 3L)
        ", ", if (n > 5L)
        " ..., ", numform(xx[i2]), "\n", sep = "")
        cat("\n")
    }
    if (label %in% c("Normal approximation",
                     "Normal Power approximation"))
        cat(attr(x, "source"), "\n")
    invisible(x)
}

plot.aggregateDist <- function(x, xlim,
                               ylab = expression(F[S](x)),
                               main = "Aggregate Claim Amount Distribution",
                               sub = comment(x), ...)
{
    ## Function plot() is used for the step cdfs and function curve()
    ## in the continuous cases.
    if ("stepfun" %in% class(x))
    {
        ## Method for class 'ecdf' will most probably be used.
        NextMethod(main = main, ylab = ylab, ...)
    }
    else
    {
        ## Limits for the x-axis are supplied if none are given
        ## in argument.
        if (missing(xlim))
        {
            mean <- get("mean", envir = environment(x))
            sd <- sqrt(get("variance", envir = environment(x)))
            xlim <- c(mean - 3 * sd, mean + 3 * sd)
        }
        curve(x, main = main, ylab = ylab, xlim = xlim, ylim = c(0, 1), ...)
    }
    mtext(sub, line = 0.5)
}

summary.aggregateDist <- function(object, ...)
    structure(object, class = c("summary.aggregateDist", class(object)), ...)

print.summary.aggregateDist <- function(x, ...)
{
    cat(ifelse(comment(x) %in%
               c("Normal approximation", "Normal Power approximation"),
               "Aggregate Claim Amount CDF:\n",
               "Aggregate Claim Amount Empirical CDF:\n"))
    q <- quantile(x, p = c(0.25, 0.5, 0.75))
    expectation <- mean(x)

    if (comment(x) %in% c("Normal approximation", "Normal Power approximation"))
    {
        min <- 0
        max <- NA
    }
    else
    {
        max <- tail(eval(expression(x), environment(x)), 1)
        min <- head(eval(expression(x), environment(x)), 1)
    }
    res <- c(min, q[c(1, 2)], expectation, q[3], max)
    names(res) <- c("Min.", "1st Qu.", "Median", "Mean", "3rd Qu.", "Max.")
    print(res, ...)
    invisible(x)
}

mean.aggregateDist <- function(x, ...)
{
    label <- comment(x)

    ## Simply return the value of the true mean given in argument in
    ## the case of the Normal and Normal Power approximations.
    if (label %in%
        c("Normal approximation", "Normal Power approximation"))
        return(get("mean", envir = environment(x)))

    ## For the recursive, exact and simulation methods, compute the
    ## mean from the stepwise cdf using the pmf saved in the
    ## environment of the object.
    drop(crossprod(get("x", envir = environment(x)),
                   get("fs", envir = environment(x))))
}

diff.aggregateDist <- function(x, ...)
{
    label <- comment(x)

    ## The 'diff' method is defined for the recursive, exact and
    ## simulation methods only.
    if (label == "Normal approximation" || label == "Normal Power approximation")
        stop("function not defined for approximating distributions")

    ## The probability vector is already stored in the environment of
    ## the "aggregateDist" object.
    get("fs", environment(x))
}
