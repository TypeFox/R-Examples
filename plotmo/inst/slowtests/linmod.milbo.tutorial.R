# linmod.milbo.tutorial.R:
#
# linmod code from Stephen Milborrow "Guidelines for S3 Regression Models"

linmod <- function(...) UseMethod("linmod")

linmod.fit <- function(x, y) # internal function, not for the casual user
{                            # first column of x is the intercept (all 1s)

    qx <- qr(x)                         # QR-decomposition of x
    y <- as.vector(as.matrix(y))        # necessary when y is a data.frame
    coef <- solve.qr(qx, y)             # compute (x'x)^(-1) x'y
    df.residual <- nrow(x) - ncol(x)    # degrees of freedom
    sigma2 <- sum((y - x %*% coef)^2) / df.residual  # variance of residuals
    vcov <- sigma2 * chol2inv(qx$qr)    # covar mat is sigma^2 * (x'x)^(-1)
    colnames(vcov) <- rownames(vcov) <- colnames(x)
    fitted.values <- qr.fitted(qx, y)

    fit <- list(coefficients  = coef,
                residuals     = y - fitted.values,
                fitted.values = fitted.values,
                vcov          = vcov,
                sigma         = sqrt(sigma2),
                df.residual   = df.residual)

    class(fit) <- "linmod"
    fit
}
linmod.default <- function(x, y, ...)
{
    fit <- linmod.fit(cbind("(Intercept)"=1, as.matrix(x)), y)
    fit$call <- match.call()
    fit
}
linmod.formula <- function(formula, data=parent.frame(), ...)
{
    mf <- model.frame(formula=formula, data=data)
    terms <- attr(mf, "terms")
    fit <- linmod.fit(model.matrix(terms, mf), model.response(mf))
    fit$terms <- terms
    fit$call <- match.call()
    fit
}
predict.linmod <- function(object, newdata=NULL, ...)
{
    if(is.null(newdata))
        y <- fitted(object)
    else {
        if(is.null(object$terms))              # x,y interface
            x <- cbind(1, as.matrix(newdata))  # columns must be in same order as orig x
        else {                                 # formula interface
            terms <- delete.response(object$terms)
            x <- model.matrix(terms, model.frame(terms, as.data.frame(newdata)))
        }
        y <- as.vector(x %*% coef(object))
    }
    y
}
