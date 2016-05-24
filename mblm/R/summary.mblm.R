"summary.mblm" <-
function (object, ...) 
{
    z <- object
    p <- z$rank
    Qr <- object$qr
    if (is.null(z$terms) || is.null(Qr)) 
        stop("invalid 'lm' object:  no 'terms' nor 'qr' component")
    n <- NROW(Qr$qr)
    rdf <- n - p
    if (is.na(z$df.residual) || rdf != z$df.residual) 
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    p1 <- 1:p
    r <- z$residuals
    f <- z$fitted
    w <- z$weights
    if (is.null(w)) {
        mss <- if (attr(z$terms, "intercept")) 
            sum((f - mean(f))^2)
        else sum(f^2)
        rss <- sum(r^2)
    }
    else {
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
    }
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    madval <- c(mad(z$intercepts),mad(z$slopes));
    est <- z$coefficients[Qr$pivot[p1]]
    vval <- c(wilcox.test(z$intercepts)$statistic,wilcox.test(z$slopes)$statistic);
    pval <- c(wilcox.test(z$intercepts)$p.value,wilcox.test(z$slopes)$p.value);
    ans <- z[c("call", "terms")]
    ans$residuals <- r;
    ans$coefficients <- cbind(est, madval, vval, pval)
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]], 
        c("Estimate", "MAD", "V value", "Pr(>|V|)"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    class(ans) <- "summary.lm"
    ans
}

