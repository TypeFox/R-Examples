# linmod.R: linmod with checks suitable for production code
#
# The original minimal version of linmod is in
# Friedrich Leisch "Creating R Packages: A Tutorial".
#
# Please see also www.milbo.org/doc/modguide.pdf.
#
# Tests for this code may be found in test.ltut.bat in the plotmo package.

linmod <- function(...) UseMethod("linmod")

# internal function, not for the casual user
# first column of x is the intercept (all 1s)

linmod.fit <- function(x = stop("no 'x' argument"),
                       y = stop("no 'y' argument"),
                       ...)
{
    stop.if.dot.arg.used(...)
    x <- check.linmod.x(x)
    y <- check.linmod.y(x, y)

    qx <- qr(x)                            # QR-decomposition of x
    if(qx$rank < ncol(x))
        stop("'x' is singular (it has ", ncol(x),
             " columns but its rank is ", qx$rank, ")")
    coef <- solve.qr(qx, y)                # compute (x'x)^(-1) x'y
    stopifnot(!anyNA(coef))                # should never happen
    df.residual <- nrow(x) - ncol(x)       # degrees of freedom
    stopifnot(df.residual > 0) # should have been caught by singular check above
    sigma2 <- sum((y - x %*% coef)^2) / df.residual # variance of residuals
    vcov <- sigma2 * chol2inv(qx$qr)       # covar mat is sigma^2 * (x'x)^(-1)
    fitted.values <- qr.fitted(qx, y)

    colnames(vcov) <- rownames(vcov) <- colnames(x)
    names(fitted.values) <- rownames(x)
    colnames(coef) <- colnames(y)

    fit <- list(coefficients  = coef,
                residuals     = y - fitted.values,
                fitted.values = fitted.values,
                vcov          = vcov,
                sigma         = sqrt(sigma2),
                df.residual   = df.residual)

    class(fit) <- "linmod"
    fit
}
check.linmod.x <- function(x)
{
    if(!is.matrix(x))
        stop("'x' is not a matrix (or could not be converted to a matrix)")
    if(NROW(x) == 0 || NCOL(x) == 0)
        stop("'x' is empty")
    if(anyNA(x))
        stop("NA in 'x'") # TODO ideally we should say where in x the NA is
    # checking just the first column of x suffices because x is a matrix
    # is.logical allowed because qr etc. know how to deal with logical vars
    if(!is.numeric(x[,1]) && !is.logical(x[,1]))
        stop("non-numeric column in 'x'")
    # the model must have an intercept (TODO until the test suite is extended)
    if(any(x[,1] != 1))
        stop("the first column of 'x' is not an intercept column (all 1s)")

    # ensure all columns in x are named (needed for names in vcov etc.)
    # use the same naming convention as lm (prefix for unnamed cols is "V")
    missing.colnames <-
        if(is.null(colnames(x))) 1:NCOL(x)
        else                     nchar(colnames(x)) == 0
    colnames(x)[missing.colnames] <-
        c("(Intercept)",
          paste("V", seq_len(NCOL(x)-1), sep=""))[missing.colnames]
    duplicated <- which(duplicated(colnames(x)))
    if(length(duplicated))
        stop("column name \"", colnames(x)[duplicated[1]],
             "\" in 'x' is duplicated")
    x
}
check.linmod.y <- function(x, y)
{
    # as.vector(as.matrix(y)) is necessary when y is a data.frame
    # (as.vector alone on a data.frame returns a data.frame)
    # we transform y in multiple steps for a clearer traceback() if error
    y <- as.matrix(y)
    y <- as.vector(y)
    if(length(y) == 0)
        stop("'y' is empty")
    if(anyNA(y))
        stop("NA in 'y'")
    if(!is.numeric(y) && !is.logical(y))
        stop("'y' is not numeric or logical")
    if(length(y) != nrow(x))
        stop("nrow(x) is ", nrow(x), " but length(y) is ", length(y))
    y
}
linmod.default <- function(x = stop("no 'x' argument"),
                           y = stop("no 'y' argument"),
                           keep = FALSE,
                           ...)
{
    stop.if.dot.arg.used(...)
    x.original <- x
    x <- as.matrix(x)
    # use name "(Intercept)" here so coef names match linmod.formula
    x <- cbind("(Intercept)"=1, x)
    fit <- linmod.fit(x, y)
    fit$call <- match.call()
    if(keep) {
        fit$x <- x.original
        fit$y <- y
    }
    fit
}
linmod.formula <- function(formula = stop("no 'formula' argument"),
                           data = parent.frame(),
                           ...)
{
    stop.if.dot.arg.used(...)
    # note that na.action=na.pass because we will catch NAs later
    # in linmod.fit, for uniformity with linmod.default
    mf <- model.frame(formula=formula, data=data, na.action=na.pass)
    terms <- attr(mf, "terms")
    x <- model.matrix(terms, mf)
    y <- model.response(mf)
    fit <- linmod.fit(x, y)
    fit$terms <- terms
    fit$call <- match.call()
    fit
}
predict.linmod <- function(object = stop("no 'object' argument"),
                           newdata = NULL,
                           ...)
{
    # following commented because by default plotmo passes a type arg to predict
    # stop.if.dot.arg.used(...)

    if(is.null(newdata))
        y <- fitted(object)
    else{
        if(NROW(newdata) == 0)
            stop("'newdata' is empty")
        if(is.null(object$terms)) {            # x,y interface
            x <- as.matrix(newdata) # columns must be in the same order as orig x
            x <- cbind(1, x)
        } else {                                # formula interface
            terms <- delete.response(object$terms)
            # TODO The following code can issue quite obscure
            #      error messages for bad newdata.  For example
            #          predict(obj, newdata=1:3)
            #      causes
            #          eval(expr, envir, enclos) : object 'varname' not found
            #
            # na.action=na.pass because we will catch NAs a little later

            x <- as.data.frame(newdata)
            x <- model.frame(terms, x, na.action=na.pass)
            x <- model.matrix(terms, x)
        }
        # TODO The following tests suffice to catch all incorrect input (I believe),
        # but aren't ideal in that they don't always direct you to the root cause
        # of the problem.  For example, strings in newdata that get converted to
        # factors by model.matrix can cause the wrong number of columns in x.
        #
        # TODO if x has colnames, warn if not the same as original colnames

        if(ncol(x) != length(object$coefficients))
            stop("ncol(newdata) is ", ncol(x)-1, " but should be ",
                 length(object$coefficients)-1) # -1 for intercept
        if(anyNA(x))
            stop("NA in 'newdata'")
        if(!is.numeric(x[,1]) && !is.logical(x[1,]))
            stop("non-numeric column in 'newdata'")
        y <- as.vector(x %*% coef(object))
    }
    y
}
summary.linmod <- function(object = stop("no 'object' argument"), ...)
{
    stop.if.dot.arg.used(...)
    se <- sqrt(diag(object$vcov))
    t.value <- coef(object) / se

    coefficients <- cbind(Estimate = coef(object),
                          StdErr   = se,
                          t.value  = t.value,
                          p.value  = 2 * pt(-abs(t.value), df=object$df))

    res <- list(call         = object$call,
                coefficients = coefficients)

    class(res) <- "summary.linmod"
    res
}
print.linmod <- function(x = stop("no 'x' argument"), ...)
{
    stop.if.dot.arg.used(...)
    cat("Call: ") # lm has a newline here, but a space is more compact
    cat(strwrap(paste0(deparse(x$call, control=NULL, nlines=5),
                       sep=" ", collapse=" "), exdent=6), sep="\n")
    cat("\n")
    print(x$coefficients)
    x
}
print.summary.linmod <- function(x = stop("no 'x' argument"), ...)
{
    stop.if.dot.arg.used(...)
    cat("Call: ")
    cat(strwrap(paste0(deparse(x$call, control=NULL, nlines=5),
                       sep=" ", collapse=" "), exdent=6), sep="\n")
    cat("\n")
    printCoefmat(x$coefficients, P.value=TRUE, has.Pvalue=TRUE)
    x
}
# stop.if.dot.arg.used will cause an error message if any args are passed to it.
# We use it to test if any dots arg of the calling function was used, for
# functions that must have a dots arg (to match the generic method) but
# don't actually use the dots.
#
# TODO R version 3.3-0 will have a function chkDots which should be used instead

stop.if.dot.arg.used <- function()
{
    NULL
}
