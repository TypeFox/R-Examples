MISE = function (actual, estimate, neval = 1000)
{
    m <- frequency(actual$time)
    s <- start(actual$time)
    x <- actual$x
    p <- length(x)
    n <- ncol(actual$y)
    if (length(estimate$x) != p | nrow(actual$y) != p | p !=
        nrow(estimate$y) | n != ncol(estimate$y))
        stop("Dimensions of inputs don't match")
    if (max(abs(actual$x - estimate$x)) > 1e-05)
        stop("Different x variables")
    e <- estimate$y - actual$y
    e.big <- pe.big <- matrix(NA, nrow = neval, ncol = n)
    calcpc <- TRUE
    if (calcpc) {
        pe <- e/actual$y
        pe.big <- matrix(NA, nrow = neval, ncol = n)
    }
    for (i in 1:n) {
        if (sum(is.na(e[, i])) == 0) {
            idx <- (abs(e[, i]) == Inf) | is.na(e[, i])
            e.big[, i] <- spline(x[!idx], e[!idx, i], n = neval,
                method = "natural")$y
            if (calcpc) {
                idx <- (abs(pe[, i]) == Inf) | is.na(pe[, i])
                pe.big[, i] <- spline(x[!idx], pe[!idx, i], n = neval,
                  method = "natural")$y
            }
        }
    }
    delta <- (max(x) - min(x))/neval
    out <- list(x = x, error = e)
    out$MIE <- ts(colSums(e.big, na.rm = TRUE) * delta, start = s,
        frequency = m)
    out$MIAE <- ts(colSums(abs(e.big), na.rm = TRUE) * delta,
        start = s, frequency = m)
    out$MISE <- ts(colSums(e.big^2, na.rm = TRUE) * delta, start = s,
        frequency = m)
    out$ME <- rowMeans(e, na.rm = TRUE)
    out$MAE <- rowMeans(abs(e), na.rm = TRUE)
    out$MSE <- rowMeans(e^2, na.rm = TRUE)
    if (calcpc) {
        out$MIPE <- ts(colSums(pe.big, na.rm = TRUE) * delta,
            start = s, frequency = m)
        out$MIAPE <- ts(colSums(abs(pe.big), na.rm = TRUE) *
            delta, start = s, frequency = m)
        out$MPE <- rowMeans(pe, na.rm = TRUE)
        out$MAPE <- rowMeans(abs(pe), na.rm = TRUE)
    }
    return(out)
}
