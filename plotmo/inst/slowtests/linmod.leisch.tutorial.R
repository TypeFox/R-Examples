# linmod.leisch.tutorial.R:
#
# linmod code from Friedrich Leisch "Creating R Packages: A Tutorial"

linmodEst <- function(x, y)
{
    ## compute QR-decomposition of x
    qx <- qr(x)

    ## compute (x'x)^(-1) x'y
    coef <- solve.qr(qx, y)

    ## degrees of freedom and standard deviation of residuals
    df <- nrow(x)-ncol(x)
    sigma2 <- sum((y - x%*%coef)^2)/df

    ## compute sigma^2 * (x'x)^-1
    vcov <- sigma2 * chol2inv(qx$qr)
    colnames(vcov) <- rownames(vcov) <- colnames(x)

    list(coefficients = coef,
         vcov = vcov,
         sigma = sqrt(sigma2),
         df = df)
}
print.linmod <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\nCoefficients:\n")
    print(x$coefficients)
}
summary.linmod <- function(object, ...)
{
    se <- sqrt(diag(object$vcov))
    tval <- coef(object) / se

    TAB <- cbind(Estimate = coef(object),
                 StdErr = se,
                 t.value = tval,
                 p.value = 2*pt(-abs(tval), df=object$df))

    res <- list(call=object$call,
                coefficients=TAB)

    class(res) <- "summary.linmod"
    res
}
print.summary.linmod <- function(x, ...)
{
    cat("Call:\n")
    print(x$call)
    cat("\n")
    printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
}

linmod <- function(x, ...) UseMethod("linmod")

linmod.default <- function(x, y, ...)
{
    x <- as.matrix(x)
    y <- as.numeric(y)

    est <- linmodEst(x, y)

    est$fitted.values <- as.vector(x %*% est$coefficients)
    est$residuals <- y - est$fitted.values
    est$call <- match.call()

    class(est) <- "linmod"
    est
}
linmod.formula <- function(formula, data=list(), ...)
{
    mf <- model.frame(formula=formula, data=data)
    x <- model.matrix(attr(mf, "terms"), data=mf)
    y <- model.response(mf)

    est <- linmod.default(x, y, ...)
    est$call <- match.call()
    est$formula <- formula
    est
}
predict.linmod <- function(object, newdata=NULL, ...)
{
    if(is.null(newdata))
        y <- fitted(object)
    else{
        if(!is.null(object$formula)){
            ## model has been fitted using formula interface
            x <- model.matrix(object$formula, newdata)
        } else{
            x <- newdata
        }
        y <- as.vector(x %*% coef(object))
    }
    y
}
