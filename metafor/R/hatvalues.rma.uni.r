hatvalues.rma.uni <-
function (model, type = "diagonal", ...) 
{
    if (!is.element("rma.uni", class(model))) 
        stop("Argument 'model' must be an object of class \"rma.uni\".")
    na.act <- getOption("na.action")
    if (!is.element(na.act, c("na.omit", "na.exclude", "na.fail", 
        "na.pass"))) 
        stop("Unknown 'na.action' specified under options().")
    type <- match.arg(type, c("diagonal", "matrix"))
    x <- model
    if (x$weighted) {
        if (is.null(x$weights)) {
            wi <- 1/(x$vi + x$tau2)
            W <- diag(wi, nrow = x$k, ncol = x$k)
            H <- x$X %*% (x$vb/x$s2w) %*% crossprod(x$X, W)
        }
        else {
            A <- diag(x$weights, nrow = x$k, ncol = x$k)
            stXAX <- .invcalc(X = x$X, W = A, k = x$k)
            H <- x$X %*% stXAX %*% crossprod(x$X, A)
        }
    }
    else {
        stXX <- .invcalc(X = x$X, W = diag(x$k), k = x$k)
        H <- x$X %*% tcrossprod(stXX, x$X)
    }
    if (type == "diagonal") {
        hii <- rep(NA_real_, x$k.f)
        hii[x$not.na] <- diag(H)
        hii[hii > 1 - 10 * .Machine$double.eps] <- 1
        names(hii) <- x$slab
        if (na.act == "na.omit") 
            hii <- hii[x$not.na]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(hii)
    }
    if (type == "matrix") {
        Hfull <- matrix(NA_real_, nrow = x$k.f, ncol = x$k.f)
        Hfull[x$not.na, x$not.na] <- H
        rownames(Hfull) <- x$slab
        colnames(Hfull) <- x$slab
        if (na.act == "na.omit") 
            Hfull <- Hfull[x$not.na, x$not.na, drop = FALSE]
        if (na.act == "na.fail" && any(!x$not.na)) 
            stop("Missing values in results.")
        return(Hfull)
    }
}
