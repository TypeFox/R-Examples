cm <- function(formula1,formula2, data, ...) UseMethod("cm")

cm.default <- function (formula1, formula2, data,...)
{
	cmEst <- function(x1, x2, y1, y2)
	{    
	    xpm <- (x1+x2)/2
	    xh <- x2-x1
	    ypm <- (y1+y2)/2
	    yh <- y2-y1
	    qxpm <- qr(xpm)
	    coef <- solve.qr(qxpm, ypm)
	    df <- nrow(xpm)-ncol(xpm)
	    sigma2 <- sum((ypm - xpm%*%coef)^2)/df
	    vcov <- sigma2 * chol2inv(qxpm$qr)
	    colnames(vcov) <- rownames(vcov) <- colnames(xpm)
	    list(coefficients = coef,
	        vcov = vcov,
	        sigma = sqrt(sigma2),
	        df = df)
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
    est <- cmEst(x1, x2, y1, y2)
    est$fitted.values.l <- as.vector(x1%*%est$coefficients)
    est$fitted.values.u <- as.vector(x2%*%est$coefficients)
    est$residuals.l <- y1 - est$fitted.values.l
    est$residuals.u <- y2 - est$fitted.values.u
    est$call <- match.call()
    class(est) <- "cm"
    est
}

print.cm <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    names(x$coefficients)[-1] <- seq(1,length(names(x$coefficients)[-1]))
    print(list(coefficients = x$coefficients,
	   sigma = x$sigma, df = x$df,
	   fitted.values.l = x$fitted.values.l, fitted.values.u = x$fitted.values.u,
	   residuals.l = x$residuals.l, residuals.u = x$residuals.u))
}

summary.cm <- function(object, ...)
{
    rmse.l <- sqrt(mean(object$residuals.l^2))
    rmse.u <- sqrt(mean(object$residuals.u^2))
    se <- sqrt(diag(object$vcov))
    tval <- coef(object) / se
    TAB <- cbind(Estimate = coef(object),
        StdErr = se)
    res <- list(call=object$call,
        coefficients=TAB,
        RMSE.l = rmse.l,
        RMSE.u = rmse.u)
    class(res) <- "summary.cm"
    res
}

print.summary.cm <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    rownames(x$coefficients)[-1] <- seq(1,length(rownames(x$coefficients)[-1]))
    print(x$coefficients)
    cat("\n")
    cat("RMSE.L:\n")
    print(x$RMSE.l)
    cat("RMSE.U:\n")
    print(x$RMSE.u)
}

fitted.cm <- function(object, ...)
{
    fit.Min <- object$fitted.values.l
    fit.Max <- object$fitted.values.u
    ftd <- cbind(fit.Min,
        fit.Max)
    fitted <- round(ftd, digits=3)
    class(fitted) <- "fitted.cm"
    fitted
}

residuals.cm <- function (object, ...)
{
    resid.Min <- object$residuals.l	
    resid.Max <- object$residuals.u
    resi <- cbind(resid.Min,resid.Max)
    resi <- round(resi,digits=3)
    class(resi) <- "residuals.cm"
    resi
}

cm.formula <- function(formula1,formula2,data=list(),...)
{
    mf1 <- model.frame(formula=formula1,data=data)
    x1 <- model.matrix(attr(mf1, "terms"), data=mf1)
    y1 <- model.response(mf1)
    mf2 <- model.frame(formula=formula2,data=data)
    x2 <- model.matrix(attr(mf2, "terms"), data=mf2)
    y2 <- model.response(mf2)
    est <- cm.default(formula1, formula2, data,...)
    est$call <- match.call()
    est$formula1 <- formula1
    est$formula2 <- formula2
    est
}
