treg <- function(lm.object, r, verbose=TRUE){
    if (class(lm.object)!="lm") stop("model must be class ``lm''")
    X <- model.matrix(lm.object)
  y <- lm.object$model$y
  nu <- length(y)
  w <- rep(1,nu)
  d <- 0
  convergence.criteria <- 0.0001
  converged <- FALSE
  fit <- lm.object
  while (!converged) {
    rss <- sum(w*resid(fit)^2)
    phi <- (r + 1) * rss/nu
    sigma <- sqrt(phi/r)
    t <- resid(fit)/sqrt(phi)
    w <- 1/(1+t^2)
    e <- -2*(nu*( lgamma((r+1)/2) -0.5*log(pi) - lgamma(r/2) -
                0.5*log(phi)) + (r + 1)/2*sum(log(w)))
    if (verbose)cat("-2 log Lmax =",e," with sigma = ",sigma," and scale parameter ",r,"\n")
    converged <- abs(d-e) < convergence.criteria
    d <- e
    fit <- lm.wfit(X,y,w)
#    cat(r,"\n")
  }
      fit <- c(lm.object,  list(weights = w, disparity=e, tcoef = coef(fit), r=r, sigma = sigma))
  class(fit) <- c("lm","treg")
  fit
  }
#
summary.treg <- function(object, ...)
#  summary.lm
#function (object, correlation = FALSE, symbolic.cor = FALSE,
#    ...)
{
    z <- object
    p <- z$rank
    Qr <- object$qr
#    if (is.null(z$terms) || is.null(Qr))
#        stop("invalid 'lm' object:  no 'terms' nor 'qr' component")
    n <- NROW(Qr$qr)
    rdf <- n - p
    if (is.na(z$df.residual) || rdf != z$df.residual)
        warning("residual degrees of freedom in object suggest this is not an \"lm\" fit")
    p1 <- 1:p
    r <- z$residuals
    f <- z$fitted
    w <- z$weights
        mss <- if (attr(z$terms, "intercept")) {
            m <- sum(w * f/sum(w))
            sum(w * (f - m)^2)
        }
        else sum(w * f^2)
        rss <- sum(w * r^2)
        r <- sqrt(w) * r
#
    resvar <- rss/rdf
    R <- chol2inv(Qr$qr[p1, p1, drop = FALSE])
    se <- sqrt(diag(R) * resvar)
    est <- z$coefficients[Qr$pivot[p1]]
    tval <- est/se
    ans <- z[c("call", "terms")]
    ans$residuals <- r
    ans$coefficients <- cbind(est, se, tval, 2 * pt(abs(tval),
        rdf, lower.tail = FALSE))
    dimnames(ans$coefficients) <- list(names(z$coefficients)[Qr$pivot[p1]],
        c("Estimate", "Std. Error", "t value", "Pr(>|t|)"))
    ans$aliased <- is.na(coef(object))
    ans$sigma <- sqrt(resvar)
    ans$df <- c(p, rdf, NCOL(Qr$qr))
    if (p != attr(z$terms, "intercept")) {
        df.int <- if (attr(z$terms, "intercept"))
            1
        else 0
        ans$r.squared <- mss/(mss + rss)
        ans$adj.r.squared <- 1 - (1 - ans$r.squared) * ((n -
            df.int)/rdf)
        ans$fstatistic <- c(value = (mss/(p - df.int))/resvar,
            numdf = p - df.int, dendf = rdf)
    }
    else ans$r.squared <- ans$adj.r.squared <- 0
    ans$cov.unscaled <- R
    dimnames(ans$cov.unscaled) <- dimnames(ans$coefficients)[c(1,
        1)]
#    if (correlation) {
#        ans$correlation <- (R * resvar)/outer(se, se)
#        dimnames(ans$correlation) <- dimnames(ans$cov.unscaled)
#        ans$symbolic.cor <- symbolic.cor
#    }
    print(ans$call)
    print(ans$coefficients)
    cat("Disparity = ",signif(z$disparity,5),"\n")
    cat("r value   = ", z$r,"\n")
    if (!is.null(z$na.action))
        ans$na.action <- z$na.action
    class(ans) <- c("summary.lm","summary.treg")
    ans
}

