ccrm <- function(formula1,formula2, data, ...) UseMethod("ccrm")

ccrm.default <- function(formula1, formula2, data, ...) {
	
	ccrmEst <- function(x.C, x.R, y.C, y.R) {
		
		qx.C <- qr(x.C)
		qx.R <- qr(x.R)
		coef.C <- solve.qr(qx.C, y.C)
	
		## calc coef.R
		Met4 <- list(y=y.R, w=rep(1,length(y.R)), X=x.R, C=matrix(0,0,0),
        		p=rep(1, ncol(x.R)), off=array(0,0), S=list(),sp=array(0,0),
         		Ain=diag(ncol(x.R)), bin=rep(0, ncol(x.R)))
		require("mgcv")
		Beta4 <- pcls(Met4)
		
		coef.R <- Beta4

		df.C <- nrow(x.C)-ncol(x.C)
		df.R <- nrow(x.R)-ncol(x.R)
		sigma2.C <- sum((y.C - x.C%*%coef.C)^2)/df.C
		sigma2.R <- sum((y.R - x.R%*%coef.R)^2)/df.R
		vcov.C <- sigma2.C * chol2inv(qx.C$qr)
		vcov.R <- sigma2.R * chol2inv(qx.R$qr)
		colnames(vcov.C) <- rownames(vcov.C) <- colnames(x.C)
		colnames(vcov.R) <- rownames(vcov.R) <- colnames(x.R)
		list(coefficients.C = coef.C,
			vcov.C = vcov.C,
			sigma.C = sqrt(sigma2.C),
			df.C = df.C,
			coefficients.R = coef.R,
			vcov.R = vcov.R,
			sigma.R = sqrt(sigma2.R),
			df.R = df.R)
	}

	## extract terms
    	mf.C <- model.frame(formula=formula1,data=data)
    	x.C <- model.matrix(attr(mf.C, "terms"), data=mf.C)
    	y.C <- model.response(mf.C)
    	mf.R <- model.frame(formula=formula2,data=data)
    	x.R <- model.matrix(attr(mf.R, "terms"), data=mf.R)
    	y.R <- model.response(mf.R)

	## calc
	x.C <- as.matrix(x.C)
    	x.R <- as.matrix(x.R)
    	y.C <- as.numeric(y.C)
    	y.R <- as.numeric(y.R)
	est <- ccrmEst(x.C, x.R, y.C, y.R)
	est$fitted.values.l <- as.vector(x.C%*%est$coefficients.C)-(as.vector(x.R%*%est$coefficients.R)/2)
	est$fitted.values.u <- as.vector(x.C%*%est$coefficients.C)+(as.vector(x.R%*%est$coefficients.R)/2)
	est$residuals.l <- (y.C - y.R/2) - est$fitted.values.l
	est$residuals.u <- (y.C + y.R/2) - est$fitted.values.u
	est$call <- match.call()
	class(est) <- "ccrm"
	est

}

print.ccrm <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    names(x$coefficients.C)[-1] <- seq(1,length(names(x$coefficients.C)[-1]))
    names(x$coefficients.R) <- seq(0,length(x$coefficients.R)-1[])
    names(x$coefficients.R)[1] <- "(Intercept)"
    print(list(coefficients.C = x$coefficients.C, coefficients.R = x$coefficients.R,
	   sigma.C = x$sigma.C, sigma.R = x$sigma.R,
	   df.C = x$df.C, df.R = x$df.R,
	   fitted.values.l = x$fitted.values.l, fitted.values.u = x$fitted.values.u,
	   residuals.l = x$residuals.l, residuals.u = x$residuals.u))
}

summary.ccrm <- function(object, ...)
{
    rmse.l <- sqrt(mean(object$residuals.l^2))
    rmse.u <- sqrt(mean(object$residuals.u^2))
    se.C <- sqrt(diag(object$vcov.C))
    se.R <- sqrt(diag(object$vcov.R))
    tval.C <- object$coefficients.C / se.C
    tval.R <- object$coefficients.R / se.R
    TAB.C <- cbind(Estimate.C = object$coefficients.C,
        StdErr.C = se.C)
    TAB.R <- cbind(Estimate.R = object$coefficients.R,
        StdErr.R = se.R)
    res <- list(call = object$call,
        coefficients.C = TAB.C,
        RMSE.l = rmse.l,
        coefficients.R = TAB.R,
        RMSE.u = rmse.u)
    class(res) <- "summary.ccrm"
    res
}

print.summary.ccrm <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    print(x$coefficients.C)
    cat("\n")
    print(x$coefficients.R)
    cat("\n")
    cat("RMSE.L:\n")
    print(x$RMSE.l)
    cat("\n")
    cat("RMSE.U:\n")
    print(x$RMSE.u)
}

coef.ccrm <- function(object, ...)
{
    coef.C <- object$coefficients.C
    coef.R <- object$coefficients.R
    coef <- list(coefficients.C = coef.C,
	   coefficients.R = coef.R)
    class(coef) <- "coef.ccrm"
    coef
}

print.coef.ccrm <- function(x, ...)
{
    names(x$coefficients.C)[-1] <- seq(1,length(names(x$coefficients.C)[-1]))
    names(x$coefficients.R) <- seq(0,length(x$coefficients.R)-1[])
    names(x$coefficients.R)[1] <- "(Intercept)"
    print(list(coefficients.C = x$coefficients.C,
		   coefficients.R = x$coefficients.R))
}

fitted.ccrm <- function(object, ...)
{
    fit.Min <- object$fitted.values.l
    fit.Max <- object$fitted.values.u
    ftd <- cbind(fit.Min,
        fit.Max)
    fitted <- round(ftd,digits=3)
    class(fitted) <- "fitted.ccrm"
    fitted
}

residuals.ccrm <- function (object, ...)
{
    resid.Min <- object$residuals.l
    resid.Max <- object$residuals.u
    resi <- cbind(resid.Min,resid.Max)
    resi = round(resi,digits=3)
    class(resi) <- "residuals.ccrm"
    resi
}

ccrm.formula <- function(formula1,formula2,data=list(),...)
{
    mf.C <- model.frame(formula=formula1,data=data)
    x.C <- model.matrix(attr(mf.C, "terms"), data=mf.C)
    y.C <- model.response(mf.C)
    mf.R <- model.frame(formula=formula2,data=data)
    x.R <- model.matrix(attr(mf.R, "terms"), data=mf.R)
    y.R <- model.response(mf.R)
    est <- ccrm.default(formula1, formula2, data,...)
    est$call <- match.call()
    est$formula1 <- formula1
    est$formula2 <- formula2
    est
}