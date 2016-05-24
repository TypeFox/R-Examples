ibs <-
function (x, df = NULL, knots = NULL, intercept = FALSE, Boundary.knots = range(x), 
                 from = 0, weight.fun = NULL, ...) {
    if (!is.null(weight.fun) && !is.function(weight.fun))
        stop("'weight.fun' must be a function.\n")
    bs.x <- if (is.null(knots)) {
        bs(x, df = df, intercept = intercept, Boundary.knots = Boundary.knots)
    } else {
        bs(x, knots = knots, intercept = intercept, Boundary.knots = Boundary.knots)
    } 
    kn <- attr(bs.x, "knots")
    Bkn <- attr(bs.x, "Boundary.knots")
    wk <- gaussKronrod(15)$wk
    sk <- gaussKronrod(15)$sk
    P1 <- (x + from) / 2
    P2 <- (x - from) / 2
    st <- outer(P2, sk) + P1
    out <- vector("list", 15)
    for (i in 1:15) {
        out[[i]] <- wk[i] * bs(st[, i], knots = kn, Boundary.knots = Bkn, intercept = intercept)
        if (!is.null(weight.fun)) {
            ww <- weight.fun(st[, i], x, ...)
            out[[i]] <- out[[i]] * ifelse(is.finite(ww), ww, 0)
        }
        
    }
    out <- P2 * Reduce("+", out)
    attr(out, "from") <- from
    attr(out, "weight.fun") <- weight.fun
    attr(out, "class") <- c("ibs", "basis", "matrix")
    out
}
