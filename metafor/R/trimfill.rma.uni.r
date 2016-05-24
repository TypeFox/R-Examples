trimfill.rma.uni <-
function (x, side, estimator = "L0", maxiter = 100, verbose = FALSE, 
    ...) 
{
    if (!is.element("rma.uni", class(x))) 
        stop("Argument 'x' must be an object of class \"rma.uni\".")
    if (!x$int.only) 
        stop("Trim-and-fill method only applicable for models without moderators.")
    if (missing(side)) 
        side <- NULL
    estimator <- match.arg(estimator, c("L0", "R0", "Q0"))
    if (x$k == 1) 
        stop("Stopped because k = 1.")
    yi <- x$yi
    vi <- x$vi
    weights <- x$weights
    ni <- x$ni
    if (is.null(side)) {
        res <- rma(yi, vi, weights = weights, mods = sqrt(vi), 
            intercept = TRUE, method = x$method, weighted = x$weighted, 
            ...)
        if (res$b[2] < 0) {
            side <- "right"
        }
        else {
            side <- "left"
        }
    }
    else {
        side <- match.arg(side, c("left", "right"))
    }
    if (side == "right") {
        yi <- -1 * yi
    }
    idix <- sort(yi, index.return = TRUE)$ix
    yi <- yi[idix]
    vi <- vi[idix]
    weights <- weights[idix]
    ni <- ni[idix]
    k <- length(yi)
    k0.sav <- -1
    k0 <- 0
    iter <- 0
    while (abs(k0 - k0.sav) > 0) {
        k0.sav <- k0
        iter <- iter + 1
        if (iter > maxiter) 
            stop("Trim and fill algorithm did not converge.")
        yi.t <- yi[1:(k - k0)]
        vi.t <- vi[1:(k - k0)]
        weights.t <- weights[1:(k - k0)]
        res <- rma(yi.t, vi.t, weights = weights.t, intercept = TRUE, 
            method = x$method, weighted = x$weighted, ...)
        b <- c(res$b)
        yi.c <- yi - b
        yi.c.r <- rank(abs(yi.c), ties.method = "first")
        yi.c.r.s <- sign(yi.c) * yi.c.r
        if (estimator == "R0") {
            k0 <- (k - max(-1 * yi.c.r.s[yi.c.r.s < 0])) - 1
            se.k0 <- sqrt(2 * max(0, k0) + 2)
        }
        if (estimator == "L0") {
            Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
            k0 <- (4 * Sr - k * (k + 1))/(2 * k - 1)
            varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                18 * k * k0 + 6 * k^2 * k0)
            se.k0 <- 4 * sqrt(varSr)/(2 * k - 1)
        }
        if (estimator == "Q0") {
            Sr <- sum(yi.c.r.s[yi.c.r.s > 0])
            k0 <- k - 1/2 - sqrt(2 * k^2 - 4 * Sr + 1/4)
            varSr <- 1/24 * (k * (k + 1) * (2 * k + 1) + 10 * 
                k0^3 + 27 * k0^2 + 17 * k0 - 18 * k * k0^2 - 
                18 * k * k0 + 6 * k^2 * k0)
            se.k0 <- 2 * sqrt(varSr)/sqrt((k - 1/2)^2 - k0 * 
                (2 * k - k0 - 1))
        }
        k0 <- max(0, k0)
        k0 <- round(k0)
        se.k0 <- max(0, se.k0)
        if (verbose) 
            cat("Iteration:", iter, "\tmissing =", k0, "\t  b =", 
                ifelse(side == "right", -1 * b, b), "\n")
    }
    if (k0 > 0) {
        if (side == "right") {
            yi.c <- -1 * (yi.c - b)
        }
        else {
            yi.c <- yi.c - b
        }
        yi.fill <- c(x$yi.f, -1 * yi.c[(k - k0 + 1):k])
        vi.fill <- c(x$vi.f, vi[(k - k0 + 1):k])
        weights.fill <- c(x$weights.f, weights[(k - k0 + 1):k])
        ni.fill <- c(x$ni.f, ni[(k - k0 + 1):k])
        attr(yi.fill, "measure") <- x$measure
        res <- rma(yi.fill, vi.fill, weights = weights.fill, 
            ni = ni.fill, intercept = TRUE, method = x$method, 
            weighted = x$weighted, ...)
        res$fill <- c(rep(FALSE, k), rep(TRUE, k0))
        res$ids <- c(x$ids, (x$k.f + 1):(x$k.f + k0))
        if (x$slab.null) {
            res$slab <- c(paste("Study", x$ids), paste("Filled", 
                seq_len(k0)))
        }
        else {
            res$slab <- c(x$slab, paste("Filled", seq_len(k0)))
        }
        res$slab.null <- FALSE
    }
    else {
        res <- x
        res$fill <- rep(FALSE, k)
    }
    res$k0 <- k0
    res$se.k0 <- se.k0
    res$side <- side
    res$k0.est <- estimator
    if (estimator == "R0") {
        m <- -1:(k0 - 1)
        res$p.k0 <- 1 - sum(choose(0 + m + 1, m + 1) * 0.5^(0 + 
            m + 2))
    }
    else {
        res$p.k0 <- NA
    }
    class(res) <- c("rma.uni.trimfill", class(res))
    return(res)
}
