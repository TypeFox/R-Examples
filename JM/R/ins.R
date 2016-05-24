ins <-
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), 
                 from = 0, weight.fun = NULL, ...) {
    if (!is.null(weight.fun) && !is.function(weight.fun))
        stop("'weight.fun' must be a function.\n")
    ns.x <- if (is.null(knots)) {
        ns(x, df = df, intercept = intercept, Boundary.knots = Boundary.knots)
    } else {
        ns(x, knots = knots, intercept = intercept, Boundary.knots = Boundary.knots)
    } 
    kn <- attr(ns.x, "knots")
    Bkn <- attr(ns.x, "Boundary.knots")
    wk <- gaussKronrod(15)$wk
    sk <- gaussKronrod(15)$sk
    P1 <- (x + from) / 2
    P2 <- (x - from) / 2
    st <- outer(P2, sk) + P1
    out <- vector("list", 15)
    for (i in 1:15) {
        out[[i]] <- wk[i] * ns(st[, i], knots = kn, Boundary.knots = Bkn, intercept = intercept)
        if (!is.null(weight.fun)) {
            ww <- weight.fun(st[, i], x, ...)
            out[[i]] <- out[[i]] * ifelse(is.finite(ww), ww, 0)
        }
    }
    out <- P2 * Reduce("+", out)
    attr(out, "from") <- from
    attr(out, "weight.fun") <- weight.fun
    attr(out, "class") <- c("ins", "basis", "matrix")
    out
}
