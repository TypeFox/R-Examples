
# $Id: parm.R 243 2008-07-22 16:33:38Z thothorn $

parm <- function(coef, vcov, df = 0) {

    if (length(coef) != nrow(vcov) ||
        length(coef) != ncol(vcov))
        stop("dimensions don't match")

    if (is.null(names(coef)))
        names(coef) <- paste("V", 1:length(coef), sep = "")
    if (is.null(colnames(vcov)))
        colnames(vcov) <- names(coef)
    if (is.null(rownames(vcov)))
        rownames(vcov) <- names(coef)

    if (!is.numeric(coef) || !is.vector(coef))
        stop(sQuote("coef"), " is not a numeric vector")

    if (!is.numeric(vcov) || !is.matrix(vcov))
        stop(sQuote("vcov"), " is not a numeric matrix")

    if (!isSymmetric(vcov, tol = sqrt(.Machine$double.eps)))
        stop(sQuote("vcov"), " is not a symmetric matrix")

    ret <- list(coef = coef, vcov = vcov, df = df)
    class(ret) <- "parm"
    ret
}

coef.parm <- function(object, ...)
    object$coef

vcov.parm <- function(object, ...)
    object$vcov

