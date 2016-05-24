plot.glmpath.cr <-
function (x, xvar = c("norm", "lambda", "step"), type = c("coefficients", 
    "aic", "bic"), xlimit = NULL, predictor = FALSE, 
    omit.zero = TRUE, breaks = FALSE, mar = NULL, eps = .Machine$double.eps, 
    main = NULL, ...) 
{
    object <- x
    ii <- object$new.A
	ii[length(ii)] <- TRUE
	summary.object<-summary(object)
    lam <- object$lambda[ii]
    xvar <- match.arg(xvar)
    type <- match.arg(type)
    coef.pred <- scale(object$b.predictor[ii, -1], FALSE, 1/object$sdx)
    coef.corr <- scale(object$b.corrector[ii, -1], FALSE, 1/object$sdx)
    xnames <- object$xnames[-1]
    if (omit.zero) {
        c1 <- drop(rep(1, nrow(coef.corr)) %*% abs(coef.corr))
        nonzero <- c1 > eps
        xnames <- xnames[nonzero]
        coef.pred <- coef.pred[, nonzero, drop = FALSE]
        coef.corr <- coef.corr[, nonzero, drop = FALSE]
    }
    m <- ncol(coef.pred)
    k <- nrow(coef.pred)
    s <- switch(xvar, norm = if (is.null(object$nopenalty.subset)) 
        apply(abs(coef.corr), 1, sum)
    else apply(abs(coef.corr[, -object$nopenalty.subset, drop = FALSE]), 
        1, sum), lambda = lam, step = as.numeric(gsub("Step ","",dimnames(summary.object)[[1]])))
    if (xvar != "lambda") {
        if (is.null(xlimit)) 
            xlimit <- max(s)
        else if (xlimit <= min(s)) 
            stop("Increase xlimit.")
        xi <- s <= xlimit
    }
    else {
        if (is.null(xlimit)) 
            xlimit <- min(s)
        else if (xlimit >= max(s)) 
            stop("Decrease xlimit.")
        xi <- s >= xlimit
    }
    k <- max(which(xi))
    xname <- switch(xvar, norm = "|beta|", lambda = "lambda", 
        step = "step")
    if (!is.null(mar)) 
        par(mar = mar)
    if (type == "aic") {
        aic <- summary.object$AIC
        plot(s[xi], aic, xlab = xname, ylab = "AIC", type = "b", 
            pch = 16, cex = 0.3, ...)
        if (is.null(main)) 
            title("AIC", line = 2.5)
        else title(main, line = 2.5)
    }
    else if (type == "bic") {
        bic <- summary.object$BIC
        plot(s[xi], bic, xlab = xname, ylab = "BIC", type = "b", 
            pch = 16, cex = 0.3, ...)
        if (is.null(main)) 
            title("BIC", line = 2.5)
        else title(main, line = 2.5)
    }
    else {
        ylab <- ifelse(object$standardize, "Standardized coefficients", 
            "Coefficients")
        matplot(s[xi], coef.corr[xi, ], xlab = xname, ..., type = "b", 
            pch = "*", ylab = ylab, lty = 1)
        if (is.null(main)) 
            title("Coefficient path", line = 2.5)
        else title(main, line = 2.5)
        abline(h = 0, lty = 3)
        axis(4, at = coef.corr[k, ], labels = xnames, cex = 0.8, 
            adj = 0, las = 1)
        if (predictor) {
            for (i in 1:m) segments(s[xi][-k], coef.corr[xi, 
                ][-k, i], s[xi][-1], coef.pred[xi, ][-1, i], 
                lty = 2, col = i)
        }
    }
    if (breaks) {
        new <- object$new.A[ii] & xi
        axis(3, at = s[new], labels = object$new.df[ii][new], 
            cex = 0.8)
        abline(v = s[new])
    }
}
