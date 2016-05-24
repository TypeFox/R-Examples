### implements either (8) or (41) of the design doc
### if argument zwz is supplied, does (8), otherwise does (41)

newpickle <- function(alphaceesigma, fixed, random, obj, y, origin, zwz,
    deriv = 0)
{
    stopifnot(inherits(obj, "aster"))
    if (missing(y)) {
        y <- obj$x
    } else {
        stopifnot(is.matrix(y))
        stopifnot(is.numeric(y))
        stopifnot(is.finite(y))
        stopifnot(dim(y) == dim(obj$x))
    }
    if (! missing(origin)) {
        stopifnot(is.matrix(origin))
        stopifnot(is.numeric(origin))
        stopifnot(is.finite(origin))
        stopifnot(dim(origin) == dim(obj$origin))
    }
    stopifnot(is.matrix(fixed))
    stopifnot(is.numeric(fixed))
    stopifnot(is.finite(fixed))
    nfix <- ncol(fixed)
    stopifnot(is.matrix(random) | is.list(random))
    if (! is.list(random))
        random <- list(random)
    for (i in seq(along = random)) {
        foo <- random[[i]]
        if (! is.matrix(foo))
            stop("random not matrix or list of matrices")
        if (! is.numeric(foo))
            stop("random not numeric matrix or list of such")
        if (! all(is.finite(foo)))
            stop("some random effects model matrix not all finite")
        if (nrow(foo) != nrow(fixed))
            stop("fixed and random effects model matrices with different nrow")
    }
    nrand <- sapply(random, ncol)
    if (! missing(zwz)) {
        stopifnot(is.matrix(zwz))
        stopifnot(is.numeric(zwz))
        stopifnot(is.finite(zwz))
        if (any(dim(zwz) != sum(nrand)))
            stop("zwz not square matrix with dimension = number of random effects")
    }
    stopifnot(length(deriv) == 1)
    stopifnot(deriv %in% c(0, 1))
    if (missing(zwz) & deriv != 0)
        stop("derivatives cannot be done unless zwz is supplied")
    stopifnot(is.vector(alphaceesigma))
    stopifnot(is.numeric(alphaceesigma))
    stopifnot(is.finite(alphaceesigma))
    if (length(alphaceesigma) != nfix + sum(nrand) + length(nrand))
        stop("alphaceesigma wrong length")

    idx <- seq(along = alphaceesigma)
    is.alpha <- idx <= nfix
    is.cee <- nfix < idx & idx <= nfix + sum(nrand)
    is.sigma <- nfix + sum(nrand) < idx
    alpha <- alphaceesigma[is.alpha]
    cee <- alphaceesigma[is.cee]
    sigma <- alphaceesigma[is.sigma]

    a <- as.vector(rep(sigma, times = nrand))
    bee <- a * cee

    modmat <- cbind(fixed, Reduce(cbind, random))
    ### note: despite documentation of the mlogl function, it actually
    ### works to have modmat a matrix rather than a 3-way array
    mout <- mlogl(c(alpha, bee), obj$pred, obj$fam, y, obj$root, modmat,
        deriv = 2, famlist = obj$famlist, origin = origin)

    idx.too <- seq(along = mout$gradient)
    is.alpha.too <- idx.too <= nfix
    is.cee.too <- nfix < idx.too

    if (missing(zwz)) {
        zwz <- mout$hessian
        zwz <- zwz[is.cee.too, ]
        zwz <- zwz[ , is.cee.too]
    }

    bigh <- zwz * outer(a, a) + diag(length(a))
    bigh.chol <- chol(bigh)
    val <- mout$value + sum(cee^2) / 2 + sum(log(diag(bigh.chol)))
    if (deriv == 0)
        return(list(value = val))

    pa <- mout$gradient[is.alpha.too]
    ### Z^T (y - mu^*)
    zymoo <- mout$gradient[is.cee.too]
    pc <- zymoo * a + cee

    bigh.inv <- chol2inv(bigh.chol)
    idx <- rep(seq(along = sigma), times = nrand)
    ps <- rep(NaN, length(sigma))
    for (k in seq(along = sigma)) {
        eek <- as.numeric(idx == k)
        ps[k] <- sum(bigh.inv * zwz * outer(a, eek)) + sum(zymoo * eek * cee)
    }

    return(list(value = val, gradient = c(pa, pc, ps)))
}

