
"sm.survival" <-
function (x, y, status, h, hv = 0.05, p = 0.5, status.code = 1,
    ...)
{
    opt <- sm.options(list(...))
    replace.na(opt, display, "line")
    replace.na(opt, ngrid, 50)
    replace.na(opt, xlab, deparse(substitute(x)))
    replace.na(opt, ylab, deparse(substitute(y)))
    replace.na(opt, eval.points, seq(min(x), max(x), length = opt$ngrid))
    eval.points <- opt$eval.points
    if (!(opt$display %in% "none" | opt$add == TRUE)) {
        plot(x, y, type = "n", xlab = opt$xlab, ylab = opt$ylab, ...)
        text(x[status == status.code], y[status == status.code], "x")
        text(x[status != status.code], y[status != status.code], "o")
    }
    n <- length(x)
    ne <- length(eval.points)
    xr <- x[order(y)]
    statusr <- status[order(y)]
    yr <- sort(y)
    w <- matrix(rep(eval.points, rep(n, ne)), ncol = n, byrow = TRUE)
    w <- w - matrix(rep(xr, ne), ncol = n, byrow = TRUE)
    w <- exp(-0.5 * (w/h)^2)
    wf <- t(apply(w, 1, rev))
    wf <- t(apply(wf, 1, cumsum))
    wf <- t(apply(wf, 1, rev))
    w <- w/wf
    st <- rep(0, n)
    st[statusr == status.code] <- 1
    w <- 1 - w * matrix(rep(st, ne), ncol = n, byrow = TRUE)
    w <- w[, st == 1]
    if (ne == 1)
        w <- matrix(w, ncol = length(w))
    yw <- yr[st == 1]
    w <- t(apply(w, 1, cumprod))
    w <- cbind(rep(1, ne), w)
    j <- -t(apply(w, 1, diff))
    J <- t(apply(j, 1, cumsum))
    wd <- J - p
    w <- exp(-0.5 * (wd/hv)^2)
    ns <- length(yw)
    s0 <- w %*% rep(1, ns)
    s1 <- (w * wd) %*% rep(1, ns)
    s2 <- (w * wd^2) %*% rep(1, ns)
    w <- w * (matrix(rep(s2, ns), ncol = ns) - wd * matrix(rep(s1,
        ns), ncol = ns))
    w <- w/(matrix(rep(s2, ns), ncol = ns) * matrix(rep(s0, ns),
        ncol = ns) - matrix(rep(s1, ns), ncol = ns)^2)
    estimate <- w %*% yw
    if (!(opt$display %in% "none"))
        lines(eval.points, estimate, lty = opt$lty)
    invisible(list(estimate = estimate, eval.points = eval.points,
        h = h, hv = hv, call = match.call()))
}
