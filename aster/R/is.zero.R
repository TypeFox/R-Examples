### implements test of (35) of the design doc

is.zero <- function(alphabeenu, fixed, random, obj, y, origin, zwz,
    tolerance = sqrt(.Machine$double.eps))
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
    if (missing(origin)) {
        origin <- obj$origin
    } else {
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
    stopifnot(is.matrix(zwz))
    stopifnot(is.numeric(zwz))
    stopifnot(is.finite(zwz))
    if (any(dim(zwz) != sum(nrand)))
        stop("zwz not square matrix with dimension = number of random effects")
    stopifnot(is.vector(alphabeenu))
    stopifnot(is.numeric(alphabeenu))
    stopifnot(is.finite(alphabeenu))
    if (length(alphabeenu) != nfix + sum(nrand) + length(nrand))
        stop("alphabeenu wrong length")

    idx <- seq(along = alphabeenu)
    is.alpha <- idx <= nfix
    is.bee <- nfix < idx & idx <= nfix + sum(nrand)
    is.nu <- nfix + sum(nrand) < idx
    alpha <- alphabeenu[is.alpha]
    bee <- alphabeenu[is.bee]
    nu <- alphabeenu[is.nu]
    dee <- rep(nu, times = nrand)

    if (all(nu > tolerance))
        return(rep(FALSE, length(nu)))
    if (any(nu < (- tolerance)))
        stop("apparently negative components of nu, impossible")
    nu[nu < tolerance] <- 0

    modmat <- cbind(fixed, Reduce(cbind, random))
    ### note: despite documentation of the mlogl function, it actually
    ### works to have modmat a matrix rather than a 3-way array
    mout <- mlogl(c(alpha, bee), obj$pred, obj$fam, y, obj$root, modmat,
        deriv = 2, famlist = obj$famlist, origin = origin)

    idx <- seq(along = mout$gradient)
    is.bee <- nfix < idx
    pb <- mout$gradient[is.bee]

    bigh <- sweep(zwz, 2, dee, "*") + diag(length(dee))
    bigh.inv <- solve(bigh)
    idx <- rep(seq(along = nu), times = nrand)
    pn <- rep(NaN, length(nu))
    for (k in seq(along = nu)) {
        eek <- as.numeric(idx == k)
        fook <- sweep(zwz, 2, eek, "*")
        pn[k] <- sum(t(bigh.inv) * fook) / 2
    }

    result <- rep(FALSE, length(nu))
    for (k in seq(along = nu))
        if (nu[k] == 0)
            result[k] <- pn[k] >= sum(pb[idx == k]^2) / 4
    return(result)
}

