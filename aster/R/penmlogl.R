### penmlogl is
###
###     - l(a + M alpha + Z D^(1 / 2) c) + (1 / 2) c^T c
###
### in the notation of equation (26) of inst/doc/re.tex

penmlogl <- function(parm, sigma, fixed, random, obj, y, origin) {
    stopifnot(inherits(obj, "aster"))
    stopifnot(is.numeric(fixed))
    stopifnot(is.finite(fixed))
    stopifnot(is.matrix(fixed))
    stopifnot(nrow(fixed) == length(obj$x))
    stopifnot(is.numeric(sigma))
    stopifnot(is.finite(sigma))
    stopifnot(is.numeric(parm))
    stopifnot(is.finite(parm))
    if (! is.list(random))
        random <- list(random)
    for (i in seq(along = random)) {
        if (! is.numeric(random[[i]]))
            stop("random or component of random (if list) not numeric")
        if (! all(is.finite(random[[i]])))
            stop("random or component of random (if list) not all finite")
        if (! is.matrix(random[[i]]))
            stop("random or component of random (if list) not matrix")
        if (nrow(random[[i]]) != length(obj$x))
            stop("nrow of random or component of random (if list) not length of response vector")
    }
    if (length(sigma) != length(random))
        stop("length(sigma) != length(random) if list or != 1 otherwise")
    if (length(parm) != ncol(fixed) + sum(sapply(random, ncol)))
        stop("length(parm) != number of columns of all model matrices, fixed and random")
    if (missing(y)) {
        y <- obj$x
    } else {
        stopifnot(is.numeric(y))
        stopifnot(is.matrix(y))
        stopifnot(identical(dim(y), dim(obj$root)))
        stopifnot(all(is.finite(as.vector(y))))
    }
    if (! missing(origin)) {
        stopifnot(is.numeric(origin))
        stopifnot(is.matrix(origin))
        stopifnot(identical(dim(origin), dim(obj$root)))
        stopifnot(all(is.finite(as.vector(origin))))
    }
    scalevec <- rep(1, ncol(fixed))
    penaltyvec <- rep(0, ncol(fixed))
    modmat <- fixed
    for (i in seq(along = random)) {
        foo <- sigma[i] * random[[i]]
        modmat <- cbind(modmat, foo)
        scalevec <- c(scalevec, rep(sigma[i], ncol(foo)))
        penaltyvec <- c(penaltyvec, rep(1, ncol(foo)))
    }
    ### note: despite documentation of the mlogl function, it actually
    ### works to have modmat a matrix rather than a 3-way array
    bar <- mlogl(parm, obj$pred, obj$fam, y, obj$root, modmat, deriv = 2,
        famlist = obj$famlist, origin = origin)
    val <- bar$value + sum(penaltyvec * parm^2) / 2
    grad <- bar$gradient + penaltyvec * parm
    hess <- bar$hessian + diag(penaltyvec)
    return(list(value = val, gradient = grad, hessian = hess,
        argument = parm, scale = scalevec, mlogl.hessian = bar$hessian,
        mlogl.gradient = bar$gradient))
}

penmlogl2 <- function(parm, alpha, sigma, fixed, random, obj, y, origin) {
    pout <- penmlogl(c(alpha, parm), sigma, fixed, random, obj, y, origin)
    nfix <- ncol(fixed)
    idx <- seq(along = pout$gradient) > nfix
    pout$gradient <- pout$gradient[idx]
    tmp <- pout$hessian[idx, , drop = FALSE]
    pout$hessian <- tmp[ , idx, drop = FALSE]
    pout$scale <- pout$scale[idx]
    return(pout)
}

