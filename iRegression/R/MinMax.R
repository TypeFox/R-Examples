MinMax <- function(formula1,formula2, data, ...) UseMethod("MinMax")

MinMax.default <- function(formula1,formula2, data, ...)
{
	MinMaxEst <- function(x1, x2, y1, y2)
	{
	    qx.l <- qr(x1)
	    qx.u <- qr(x2)
	    coef.l <- solve.qr(qx.l, y1)
	    coef.u <- solve.qr(qx.u, y2)
	    df.l <- nrow(x1)-ncol(x1)
	    df.u <- nrow(x2)-ncol(x2)
	    sigma2.l <- sum((y1 - x1%*%coef.l)^2)/df.l
	    sigma2.u <- sum((y2 - x2%*%coef.u)^2)/df.u
	    vcov.l <- sigma2.l * chol2inv(qx.l$qr)
	    vcov.u <- sigma2.u * chol2inv(qx.u$qr)
	    colnames(vcov.l) <- rownames(vcov.l) <- colnames(x1)
	    colnames(vcov.u) <- rownames(vcov.u) <- colnames(x2)
	    list(coefficients.l = coef.l,
	        vcov.l = vcov.l,
	        sigma.l = sqrt(sigma2.l),
	        df.l = df.l,
	        coefficients.u = coef.u,
	        vcov.u = vcov.u,
	        sigma.u = sqrt(sigma2.u),
	        df.u = df.u)
	}
    ## extract terms
    mf1 <- model.frame(formula=formula1,data=data)
    x1 <- model.matrix(attr(mf1, "terms"), data=mf1)
    y1 <- model.response(mf1)
    mf2 <- model.frame(formula=formula2,data=data)
    x2 <- model.matrix(attr(mf2, "terms"), data=mf2)
    y2 <- model.response(mf2)
    ## calc
    x1 <- as.matrix(x1)
    x2 <- as.matrix(x2)
    y1 <- as.numeric(y1)
    y2 <- as.numeric(y2)
    est <- MinMaxEst(x1, x2, y1, y2)
    est$fitted.values.l <- as.vector(x1%*%est$coefficients.l)
    est$fitted.values.u <- as.vector(x2%*%est$coefficients.u)
    est$residuals.l <- y1 - est$fitted.values.l
    est$residuals.u <- y2 - est$fitted.values.u
    est$call <- match.call()
    class(est) <- "MinMax"
    est
}

print.MinMax <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    print(list(coefficients.l = x$coefficients.l, coefficients.u = x$coefficients.u,
	   sigma.l = x$sigma.l, sigma.u = x$sigma.u,
	   df.l = x$df.l, df.u = x$df.u,
	   fitted.values.l = x$fitted.values.l, fitted.values.u = x$fitted.values.u,
	   residuals.l = x$residuals.l, residuals.u = x$residuals.u))
}

summary.MinMax <- function(object, ...)
{
    rmse.l <- sqrt(mean(object$residuals.l^2))
    rmse.u <- sqrt(mean(object$residuals.u^2))
    se.l <- sqrt(diag(object$vcov.l))
    se.u <- sqrt(diag(object$vcov.u))
    tval.l <- object$coefficients.l / se.l
    tval.u <- object$coefficients.u / se.u
    TAB.l <- cbind(Estimate.L = object$coefficients.l,
        StdErr.L = se.l)
    TAB.u <- cbind(Estimate.U = object$coefficients.u,
        StdErr.U = se.u)
    res <- list(call = object$call,
        coefficients.l = TAB.l,
        RMSE.l = rmse.l,
        coefficients.u = TAB.u,
        RMSE.u = rmse.u)
    class(res) <- "summary.MinMax"
    res
}

print.summary.MinMax <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    print(x$coefficients.l)
    cat("\n")
    cat("RMSE.L:\n")
    print(x$RMSE.l)
    cat("\n")
    print(x$coefficients.u)
    cat("\n")
    cat("RMSE.U:\n")
    print(x$RMSE.u)
}

coef.MinMax <- function(object, ...)
{
    coef.l <- object$coefficients.l
    coef.u <- object$coefficients.u
    coef <- list(coefficients.l = coef.l,
	   coefficients.u = coef.u)
    class(coef) <- "coef.MinMax"
    coef
}

print.coef.MinMax <- function(x, ...)
{
    print(list(coefficients.l = x$coefficients.l,
		   coefficients.u = x$coefficients.u))
}

fitted.MinMax <- function(object, ...)
{
    fit.Min <- object$fitted.values.l
    fit.Max <- object$fitted.values.u
    ftd <- cbind(fit.Min,
        fit.Max)
    fitted <- round(ftd,digits=3)
    class(fitted) <- "fitted.MinMax"
    fitted
}

residuals.MinMax <- function (object, ...)
{
    resid.Min <- object$residuals.l	
    resid.Max <- object$residuals.u
    resi <- cbind(resid.Min,resid.Max)
    resi <- round(resi,digits=3)
    class(resi) <- "residuals.MinMax"
    resi
}

MinMax.formula <- function(formula1,formula2,data=list(),...)
{
    mf1 <- model.frame(formula=formula1,data=data)
    x1 <- model.matrix(attr(mf1, "terms"), data=mf1)
    y1 <- model.response(mf1)
    mf2 <- model.frame(formula=formula2,data=data)
    x2 <- model.matrix(attr(mf2, "terms"), data=mf2)
    y2 <- model.response(mf2)
    est <- MinMax.default(formula1,formula2, data, ...)
    est$call <- match.call()
    est$formula1 <- formula1
    est$formula2 <- formula2
    est
}