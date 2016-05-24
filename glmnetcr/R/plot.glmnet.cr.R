plot.glmnet.cr <-
function(x, xvar = c("norm", "lambda", "step"), type = c("coefficients", 
    "aic", "bic"), omit.zero = TRUE, breaks = TRUE, mar = NULL, eps = .Machine$double.eps, 
    main = NULL, ...)
{
    object <- x
    lam <- object$lambda
    xvar <- match.arg(xvar)
    type <- match.arg(type)
    sdx<-apply(object$x,2,sd)
	beta<-as.matrix(object$beta)
    coef.pred <- scale(t(beta[1:(dim(beta)[1]-(length(levels(object$y))-1)),]), FALSE, 1/sdx)
    xnames <- dimnames(object$beta)[[1]][-c((dim(object$x)[2]+1):(dim(object$x)[2]+length(levels(object$y))-1))]
    if (omit.zero) {
        c1 <- drop(rep(1, nrow(coef.pred)) %*% abs(coef.pred))
        nonzero <- c1 > eps
        xnames <- xnames[nonzero]
        coef.pred <- coef.pred[, nonzero, drop = FALSE]
    }
    xname <- switch(xvar, norm = "|beta|", lambda = "lambda", 
        step = "step")
    s <- switch(xvar, norm = apply(abs(coef.pred), 1, sum), lambda = lam, step = seq(1,length(object$lambda)))
    phat<-predict(object)
    if (type == "aic") {
        aic <- phat$AIC
        plot(s, aic, xlab = xname, ylab = "AIC", type = "b", 
            pch = 16, cex = 0.3, ...)
        if (is.null(main)) 
            title("AIC", line = 2.5)
        else title(main, line = 2.5)
    }
    else if (type == "bic") {
        bic <- phat$BIC
        plot(s, bic, xlab = xname, ylab = "BIC", type = "b", 
            pch = 16, cex = 0.3, ...)
        if (is.null(main)) 
            title("BIC", line = 2.5)
        else title(main, line = 2.5)
    }
    else {
        ylab <-   "Coefficients"
        matplot(s, coef.pred, xlab = xname, ..., type = "b", 
            pch = "*", ylab = ylab, lty = 1)
        if (is.null(main)) 
            title("Coefficient path", line = 2.5)
        else title(main, line = 2.5)
        abline(h = 0, lty = 3)
        axis(4, at = coef.pred[length(object$lambda), ], labels = xnames, cex = 0.8, 
            adj = 0, las = 1)
     }
}

