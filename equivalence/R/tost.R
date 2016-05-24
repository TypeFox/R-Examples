tost <-
    function (x, y = NULL, epsilon = 1, paired = FALSE, var.equal = FALSE,
              conf.level = 0.95, alpha = NULL,
    ...)
{
    if (!missing(conf.level) && (length(conf.level) != 1 || !is.finite(conf.level) ||
        conf.level < 0 || conf.level > 1))
        stop("'conf.level' must be a single number between 0 and 1")
    if (!is.null(alpha)) conf.level <- 1 - alpha
    if (!is.null(y)) {
        dname <- paste(deparse(substitute(x)), "and", deparse(substitute(y)))
        if (paired)
            xok <- yok <- complete.cases(x, y)
        else {
            yok <- !is.na(y)
            xok <- !is.na(x)
        }
        y <- y[yok]
    }
    else {
        dname <- deparse(substitute(x))
        if (paired)
            stop("'y' is missing for paired test")
        xok <- !is.na(x)
        yok <- NULL
    }
    x <- x[xok]
    if (paired) {
        x <- x - y
        y <- NULL
    }
    nx <- length(x)
    mx <- mean(x)
    vx <- var(x)
    mu <- 0
    if (is.null(y)) {
        if (nx < 2)
            stop("not enough 'x' observations")
        df <- nx - 1
        stderr <- sqrt(vx/nx)
        if (stderr < 10 * .Machine$double.eps * abs(mx))
            stop("data are essentially constant")
        tstat <- (mx - mu)/stderr
        method <- if (paired)
            "Paired TOST"
        else "One Sample TOST"
        estimate <- setNames(mx, if (paired)
            "mean of the differences"
        else "mean of x")
    }
    else {
        ny <- length(y)
        if (nx < 1 || (!var.equal && nx < 2))
            stop("not enough 'x' observations")
        if (ny < 1 || (!var.equal && ny < 2))
            stop("not enough 'y' observations")
        if (var.equal && nx + ny < 3)
            stop("not enough observations")
        my <- mean(y)
        vy <- var(y)
        method <- paste(if (!var.equal)
            "Welch", "Two Sample TOST")
        estimate <- c(mx, my)
        names(estimate) <- c("mean of x", "mean of y")
        if (var.equal) {
            df <- nx + ny - 2
            v <- 0
            if (nx > 1)
                v <- v + (nx - 1) * vx
            if (ny > 1)
                v <- v + (ny - 1) * vy
            v <- v/df
            stderr <- sqrt(v * (1/nx + 1/ny))
        }
        else {
            stderrx <- sqrt(vx/nx)
            stderry <- sqrt(vy/ny)
            stderr <- sqrt(stderrx^2 + stderry^2)
            df <- stderr^4/(stderrx^4/(nx - 1) + stderry^4/(ny -
                1))
        }
        if (stderr < 10 * .Machine$double.eps * max(abs(mx),
            abs(my)))
            stop("data are essentially constant")
        tstat <- (mx - my - mu)/stderr
    }
### Implement TOST here
    pval.l <- pt(tstat, df)
    pval.u <- pt(tstat, df, lower.tail = FALSE)
    cint <- c(tstat - qt(conf.level, df),
              tstat + qt(conf.level, df)) * stderr
###
    names(tstat) <- "t"
    names(df) <- "df"
    names(mu) <- if (paired || !is.null(y))
        "difference in means"
    else "mean"
    attr(cint, "conf.level") <- conf.level
    diff <- if (!is.null(y))
        mx - my
    else estimate
    pv <- as.numeric(pt((epsilon - abs(diff)) / stderr,
                        df, lower.tail = FALSE))
    if (cint[2] < epsilon & cint[1] > -epsilon)
        result <- "rejected"
    else result <- "not rejected"
    rval <- list(parameter = df, tost.p.value = pv, result = result,
                 # conf.int = NULL,
                 tost.interval = cint, estimate = estimate, null.value = mu,
        method = method, data.name = dname, epsilon = epsilon, se.diff = stderr)
    class(rval) <- c("tost","htest")
    return(rval)
}

print.tost <- function(x, ... ) {
    NextMethod()
    cat("Epsilon:", x$epsilon, "\n")
    cat(format(100 * attr(x$tost.interval, "conf.level")),
        " percent two one-sided confidence interval (TOST interval):\n",
        " ", paste(format(c(x$tost.interval[1L], x$tost.interval[2L])),
                   collapse = " "), "\n", sep = "")
    cat("Null hypothesis of statistical difference is:", x$result, "\n")
    cat("TOST p-value:", x$tost.p.value, "\n")
    invisible(x)
}


