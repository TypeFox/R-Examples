## generate periodic B-spline basis

pbs = function (x, df = NULL, knots = NULL, degree = 3, intercept = FALSE, Boundary.knots = range(x), periodic = TRUE) {
    nx <- names(x)
    x <- as.vector(x)
    nax <- is.na(x)
    if (nas <- any(nax)) 
        x <- x[!nax]
    if (!missing(Boundary.knots)) {
        Boundary.knots <- sort(Boundary.knots)
        outside <- (ol <- x < Boundary.knots[1L]) | (or <- x > 
            Boundary.knots[2L])
    }
    else outside <- FALSE
    ord <- 1 + (degree <- as.integer(degree))
    if (ord <= 1) 
        stop("'degree' must be integer >= 1")

    if (periodic == TRUE){
        if (!missing(df) && missing(knots)) {
            nIknots <- df - 1 + (1 - intercept)
            if (nIknots < 0) {
                nIknots <- 0
                warning("'df' was too small; have used  ", ord - 
                    (1 - intercept))
            }
            knots <- if (nIknots > 0) {
                knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                    2)[-c(1, nIknots + 2)]
                stats::quantile(x[!outside], knots)
            }
        }
    } else {
        if (!missing(df) && missing(knots)) {
            nIknots <- df - ord + (1 - intercept)
            if (nIknots < 0) {
                nIknots <- 0
                warning("'df' was too small; have used  ", ord - 
                    (1 - intercept))
            }
            knots <- if (nIknots > 0) {
                knots <- seq.int(from = 0, to = 1, length.out = nIknots + 
                    2)[-c(1, nIknots + 2)]
                stats::quantile(x[!outside], knots)
            }
        }
    }
    ## get periodic B-spline knots
    getPeriodBsplineKnots = function(knots, degree=3){
        knots = sort(knots)
        nKnots = length(knots)
        if (degree < 1){
            stop("'degree' must be integer >= 1")
        }

        if (nKnots - 2 <  degree) {
            stop("number of internal knots(no boudary) must be greater than or equal to 2+ degree")
        }
        closeKnots = knots
        for (i in 1:degree){
            closeKnots = c(closeKnots, knots[nKnots] + knots[i+1] - knots[1])
        }
    
        for (i in 1:degree){
            closeKnots = c(knots[1] - knots[nKnots] + knots[nKnots-i], closeKnots)
        }
        return(closeKnots)
    }
    
    if (periodic){
        Aknots <- getPeriodBsplineKnots(sort(c(Boundary.knots, knots)), degree=degree ) ## no repeat boundary knots 3 more time
        if (any(outside)) {
            stop("some 'x' values beyond boundary knots may cause ill-conditioned bases")
        }

        basisInterior <- spline.des(Aknots, x, ord)$design
        basisInteriorLeft = basisInterior[,1:degree, drop = FALSE]
        basisInteriorRight = basisInterior[,(ncol(basisInterior)-degree+1):ncol(basisInterior), drop = FALSE]
        #print(basisInteriorLeft + basisInteriorRight)
        basis = cbind(basisInterior[,-c(1:degree, (ncol(basisInterior)-degree+1):ncol(basisInterior)), drop = FALSE],  basisInteriorLeft + basisInteriorRight)
        #df = df - degree
    } else {
        Aknots <- sort(c(rep(Boundary.knots, ord), knots))
        if (any(outside)) {
            warning("some 'x' values beyond boundary knots may cause ill-conditioned bases")
            derivs <- 0:degree
            scalef <- gamma(1L:ord)
            basis <- array(0, c(length(x), length(Aknots) - degree - 
                1L))
            if (any(ol)) {
                k.pivot <- Boundary.knots[1L]
                xl <- cbind(1, outer(x[ol] - k.pivot, 1L:degree, 
                    "^"))
                tt <- spline.des(Aknots, rep(k.pivot, ord), ord, 
                    derivs)$design
                basis[ol, ] <- xl %*% (tt/scalef)
            }
            if (any(or)) {
                k.pivot <- Boundary.knots[2L]
                xr <- cbind(1, outer(x[or] - k.pivot, 1L:degree, 
                    "^"))
                tt <- spline.des(Aknots, rep(k.pivot, ord), ord, 
                    derivs)$design
                basis[or, ] <- xr %*% (tt/scalef)
            }
            if (any(inside <- !outside)) 
                basis[inside, ] <- spline.des(Aknots, x[inside], 
                    ord)$design
        }
        else basis <- spline.des(Aknots, x, ord)$design

    }
    if (!intercept) 
        basis <- basis[, -1L, drop = FALSE]
    n.col <- ncol(basis)
    if (nas) {
        nmat <- matrix(NA, length(nax), n.col)
        nmat[!nax, ] <- basis
        basis <- nmat
    }
    dimnames(basis) <- list(nx, 1L:n.col)
    a <- list(degree = degree, knots = if (is.null(knots)) numeric(0L) else knots, 
        Boundary.knots = Boundary.knots, intercept = intercept, periodic = periodic)
    attributes(basis) <- c(attributes(basis), a)
    class(basis) <- c("pbs", "basis")
    basis

}


predict.pbs <- function(object, newx, ...)
{
    if(missing(newx))
        return(object)
    a <- c(list(x = newx), attributes(object)[
                c("degree", "knots", "Boundary.knots", "intercept", "periodic")])
    do.call("pbs", a)
}



makepredictcall.pbs <- function(var, call)
{
    if(as.character(call)[1L] != "pbs") return(call)
    at <- attributes(var)[c("degree", "knots", "Boundary.knots", "intercept", "periodic")]
    xxx <- call[1L:2]
    xxx[names(at)] <- at
    xxx
}
