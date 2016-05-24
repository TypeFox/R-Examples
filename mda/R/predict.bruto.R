predict.bruto <-
function (object, newdata, type = c("fitted", "terms"), ...)
{
    if (missing(newdata)) {
        z <- fitted(object)
        if (is.null(z)) 
            stop("need to supply newdata")
        else return(z)
    }
    d <- as.integer(dim(newdata))
    type <- match.arg(type)
    nq <- d[2]
    n <- d[1]
    if (nq != length(object$df)) 
        stop(paste("newdata should have the same number of columns",
                   "as the df component of object"))
    ybar <- object$ybar
    np <- as.integer(length(ybar))
    eta <- matrix(double(n * np), n, np)
    Type <- as.numeric(object$type)
    storage.mode(Type) <- "integer"
    storage.mode(newdata) <- "double"
    if (type == "fitted") {
        .Fortran("pbruto", newdata, n, nq, ybar, np, object$knot, object$nkmax, 
            object$nk, object$coef, Type, object$xrange, eta = eta, 
            eta, PACKAGE = "mda")$eta
    }
    else {
        ob <- as.list(seq(nq))
        names(ob) <- dimnames(newdata)[[2]]
        knot <- object$knot
        nk <- object$nk
        xrange <- object$xrange
        coef <- object$coef
        fitm <- matrix(double(n * np), n, np)
        dimnames(fitm) <- list(dimnames(newdata)[[1]], names(ybar))
        for (i in seq(nq)) {
            if (Type[i] > 1) 
                fit <- .Fortran("psspl2", newdata[, i], n, np, knot[, 
                  i], nk[i], xrange[, i], coef[, i], coef[, i], 
                  fit = fitm, as.integer(0), Type[i], PACKAGE = "mda")$fit
            else fit <- fitm
            ob[[i]] <- list(x = newdata[, i], y = fit)
        }
        ob
    }
}

