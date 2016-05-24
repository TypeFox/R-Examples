### ===== actuar: An R Package for Actuarial Science =====
###
### Minimum distance estimation.
###
### AUTHORS:  Mathieu Pigeon, Vincent Goulet <vincent.goulet@act.ulaval.ca>

mde <- function(x, fun, start, measure = c("CvM", "chi-square", "LAS"),
                weights = NULL, ...)
{
    ## General form of the function to minimize.
    myfn <- function(parm, x, weights, ...)
    {
        y <- G(parm, x, ...) - Gn(x)
        drop(crossprod(weights * y, y))
    }

    ## Extract call; used to build the call to optim().
    Call <- match.call(expand.dots = TRUE)

    ## Argument checking
    if (missing(start) || !is.list(start))
        stop("'start' must be a named list")
    if (missing(fun) || !(is.function(fun)))
        stop("'fun' must be supplied as a function")
    grouped <- inherits(x, "grouped.data")
    if (!(is.numeric(x) || grouped))
        stop("'x' must be a numeric vector or an object of class \"grouped.data\"")

    ## Make sure that any argument of 'fun' specified in '...' is held
    ## fixed.
    dots <- names(list(...))
    dots <- dots[!is.element(dots, c("upper", "lower"))]
    start <- start[!is.element(names(start), dots)]

    ## Adapt 'fun' to our needs; taken from MASS::fitdistr.
    nm <- names(start)
    f <- formals(fun)
    args <- names(f)
    m <- match(nm, args)
    if (any(is.na(m)))
        stop("'start' specifies names which are not arguments to 'fun'")
    formals(fun) <- c(f[c(1, m)], f[-c(1, m)]) # reorder arguments
    fn <- function(parm, x, ...) fun(x, parm, ...)
    if ((l <- length(nm)) > 1)
        body(fn) <- parse(text = paste("fun(x,", paste("parm[", 1:l, "]", collapse = ", "), ")"))

    measure <- match.arg(measure)

    ## Cramer-von Mises. Use the true and empirical cdf for individual
    ## data, or the true cdf and the ogive for grouped data.
    if (measure == "CvM")
    {
        G <- fn
        Gn <- if (grouped) ogive(x) else ecdf(x)
        if (is.null(weights)) weights <- 1
        Call$x <- knots(Gn)
        Call$par <- start
    }

    ## Modified Chi-square.
    if (measure == "chi-square")
    {
        if (!grouped)
            stop("\"chi-square\" measure requires an object of class \"grouped.data\"")
        if (any((nj <- x[, 2]) == 0))
            stop("frequency must be larger than 0 in all groups")
        og <- ogive(x)
        x <- knots(og)
        n <- sum(nj)
        G <- function(...) n * diff(fn(...))
        Gn <- function(...) n * diff(og(...))
        if (is.null(weights)) weights <- 1/nj
        Call$x <- x
        Call$par <- start
    }

    ## Layer average severity.
    if (measure == "LAS")
    {
        if (!grouped)
            stop("\"LAS\" measure requires an object of class \"grouped.data\"")
        e <- elev(x)
        x <- knots(e)
        G <- function(...) diff(fn(...))
        Gn <- function(...) diff(e(...))
        if (is.null(weights)) weights <- 1
        Call$x <- x
        Call$par <- start
    }

    ## optim() call
    Call[[1]] <- as.name("optim")
    Call$fun <- Call$start <- Call$measure <- NULL
    Call$fn <- myfn
    Call$weights <- weights
    Call$hessian <- FALSE
    if (is.null(Call$method))
    {
        if (any(c("lower", "upper") %in% names(Call)))
            Call$method <- "L-BFGS-B"
        else if (length(start) > 1)
            Call$method <- "BFGS"
        else
            Call$method <- "Nelder-Mead"
    }
    res <- eval(Call)

    ## Return result
    if (res$convergence > 0)
        stop("optimization failed")
    structure(list(estimate = res$par, distance = res$value),
              class = c("mde","list"))
}

print.mde <- function(x, digits = getOption("digits"), ...)
{
    ans1 <- format(x$estimate, digits = digits)
    ans1 <- sapply(ans1, function(x) paste("", x))
    nm1 <- names(ans1)
    nm1 <- paste(substring("      ", 1L, (nchar(ans1) - nchar(nm1)) %/% 2),
                 nm1)
    nm1 <- paste(nm1,
                 substring("      ", 1L, (nchar(ans1) - nchar(nm1)) %/% 2 + 1))
    names(ans1) <- nm1

    ans2 <- format(x$distance, digits = digits)
    ans2 <- sapply(ans2, function(x) paste("", x))
    nm2 <- "distance"
    nm2 <- paste(substring("      ", 1L, (nchar(ans2) - nchar(nm2)) %/% 2),
                 nm2)
    nm2 <- paste(nm2,
                 substring("      ", 1L, (nchar(ans2) - nchar(nm2)) %/% 2))
    names(ans2) <- nm2

    print(ans1, quote = FALSE)
    cat("\n")
    print(ans2, quote = FALSE)
    invisible(x)
}
