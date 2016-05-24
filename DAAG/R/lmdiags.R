lmdiags <-
function (x, which = c(1L:3L, 5L), cook.levels = c(0.5, 1), hii=NULL)
{
    dropInf <- function(x, h) {
        if (any(isInf <- h >= 1)) {
            x[isInf] <- NaN
        }
        x
    }
    if(is.null(hii))hii<- lm.influence(x, do.coef = FALSE)$hat
    out <- list(yh = NULL, rs = NULL, yhn0 = NULL, cook = NULL,
                rsp = NULL, facval=NULL)
    if (!inherits(x, "lm"))
        stop("use only with \"lm\" objects")
    if (!is.numeric(which) || any(which < 1) || any(which > 6))
        stop("'which' must be in 1:6")
    isGlm <- inherits(x, "glm")
    show <- rep(FALSE, 6)
    show[which] <- TRUE
    r <- residuals(x)
    yh <- predict(x)
    w <- weights(x)
    if (!is.null(w)) {
        wind <- w != 0
        r <- r[wind]
        yh <- yh[wind]
        w <- w[wind]
    }
    n <- length(r)
    if (any(show[2L:6L])) {
        s <- if (inherits(x, "rlm"))
            x$s
        else if (isGlm)
            sqrt(summary(x)$dispersion)
        else sqrt(deviance(x)/df.residual(x))
        if (any(show[4L:6L])) {
            cook <- if (isGlm)
                cooks.distance(x)
            else cooks.distance(x, sd = s, res = r)
        }
    }
    if (any(show[2L:3L])) {
        r.w <- if (is.null(w))
            r
        else sqrt(w) * r
        rs <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
    }
    if (any(show[5L:6L])) {
        r.hat <- range(hii, na.rm = TRUE)
        isConst.hat <- all(r.hat == 0) || diff(r.hat) < 1e-10 *
            mean(hii, na.rm = TRUE)
    }
    if (show[1L]) {
        out$yh <- yh
        out$r <- r
    }
    if (show[2L]) {
        out$rs <- rs
    }
    if (show[3L]) {
        sqrtabsr <- sqrt(abs(rs))
        yhn0 <- if (is.null(w))
            yh
        else yh[w != 0]
        if (is.null(out$rs))
            out$rs <- rs
        out$yhn0 <- yhn0
    }
    if (show[4L]) {
        out$cook <- cook
    }
    if (show[5L]) {
        r.w <- residuals(x, "pearson")
        if (!is.null(w))
            r.w <- r.w[wind]
        rsp <- dropInf(r.w/(s * sqrt(1 - hii)), hii)
        out$rsp <- rsp
    }
    if (show[6L]) {
        if (is.null(out$cook))
            out$cook <- cook
    }
    out
}
