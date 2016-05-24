fuzzy.signrank.test <- function(x,
    alternative = c("two.sided", "less", "greater"),
    mu = 0, tol = sqrt(.Machine$double.eps), alpha)
{
    alternative <- match.arg(alternative)

    if (! is.numeric(x))
        stop("'x' must be numeric")
    if (! all(is.finite(x)))
        stop("'x' must be all finite")

    if (! is.numeric(mu))
        stop("'mu' must be numeric")
    if (length(mu) != 1)
        stop("'mu' must be a single number")
    if (! is.finite(mu))
        stop("'mu' must be finite")

    if (! is.numeric(tol))
        stop("'tol' must be numeric")
    if (length(tol) != 1)
        stop("'tol' must be a single number")
    if (! is.finite(tol))
        stop("'tol' must be finite")
    if (tol < 0.0)
        stop("'tol' must be nonnegative")

    if (! missing(alpha)) {
        if (! is.numeric(alpha))
            stop("'alpha' must be numeric")
        if (length(alpha) != 1)
            stop("'alpha' must be a single number")
        if (! (is.finite(alpha) & 0 <= alpha & alpha <= 1))
            stop("'alpha' must satisfy 0 <= alpha <= 1")
    }

    dname <- deparse(substitute(x))

    sortx <- as.double(sort(x))
    out <- .Call("fpvsignrank", sortx, as.double(mu),
        alternative, as.double(tol), PACKAGE = "fuzzyRankTests")

    method <- "Wilcoxon signed rank test"
    if (missing(alpha)) {
        foo <- list(knots = out$knots, values = out$values,
            null.value = c(mu = mu), alternative = alternative,
            method = method, data.name = dname)
        return(structure(foo, class = "fuzzyranktest"))
    }

    whyknots <- out$knots
    if (alpha >= max(whyknots)) {
        reject <- 1.0
    } else if (alpha <= min(whyknots)) {
        reject <- 0.0
    } else {
        iup <- min(seq(along = whyknots)[alpha < whyknots])
        idn <- max(seq(along = whyknots)[alpha > whyknots])
        tup <- whyknots[iup] - alpha
        tdn <- alpha - whyknots[idn]
        reject <- (tdn * out$values[iup] + tup * out$values[idn]) / (tup + tdn)
    }

    foo <- list(knots = out$knots, values = out$values,
        reject.prob = reject, alpha = alpha,
        null.value = c(mu = mu), alternative = alternative,
        method = method, data.name = dname, tol = tol)
    return(structure(foo, class = "fuzzyranktest"))
}
