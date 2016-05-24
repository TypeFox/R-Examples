### minus the PQL criterion for estimating the variance components is
###
###     pout$value + log(det(I + H)) / 2
###
### where pout is the object returned by the penmlogl function, and where
###
###     H = D^(1 / 2) t(Z) W Z D^(1 / 2)
###
### (see equation (27) of inst/doc/re.tex) and I is the identity matrix.
### H is the random effects block of
###
###     pout$mlogl.hessian
###
### If lambda is the vector of eigenvalues of H, then sum(log1p(lambda))
### calculates log(det(H + I))

checkargs <- function(sigma, parm, fixed, random, obj, y, origin, cache, ...)
{
    stopifnot(inherits(obj, "aster"))
    if (! missing(y)) {
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
    if (! missing(cache))
        stopifnot(is.environment(cache))
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
        if (nrow(foo) != nrow(fixed))
            stop("fixed and random effects model matrices with different nrow")
        if (! all(is.finite(foo)))
            stop("some random effects model matrix not all finite")
    }
    nrand <- sapply(random, ncol)
    stopifnot(is.numeric(sigma))
    stopifnot(is.finite(sigma))
    if (length(sigma) != length(nrand))
        stop("length(sigma) != number of random effects model matrices")
    if ((! missing(cache)) && exists("parm", envir = cache, inherits = FALSE)) {
        stopifnot(is.numeric(cache$parm))
        stopifnot(is.finite(cache$parm))
        if (length(cache$parm) != nfix + sum(nrand))
            stop("length(cache$parm) != total number of effects (fixed and random)")
    } else {
        stopifnot(is.numeric(parm))
        stopifnot(is.finite(parm))
        if (length(parm) != nfix + sum(nrand))
            stop("length(parm) != total number of effects (fixed and random)")
    }
}

pickle <- function(sigma, parm, fixed, random, obj, y, origin, cache, ...)
{
    checkargs(sigma, parm, fixed, random, obj, y, origin, cache, ...)
    if (missing(y)) y <- obj$x
    if ((! missing(cache)) && exists("parm", envir = cache, inherits = FALSE))
        parm <- cache$parm
    trustargs <- list(objfun = penmlogl, parinit = parm, sigma = sigma,
        fixed = fixed, random = random, obj = obj, y = y, ...)
    if (is.null(trustargs$rinit)) trustargs$rinit <- 1
    if (is.null(trustargs$rmax)) trustargs$rmax <- 10
    if (! missing(origin)) trustargs$origin <- origin
    trustargs$cache <- NULL
    tout <- do.call(trust, trustargs)
    stopifnot(tout$converged)
    if ((! missing(cache))) cache$parm <- tout$argument
    nfix <- ncol(fixed)
    bigh <- tout$hessian
    bigh <- bigh[1:nrow(bigh) > nfix, , drop = FALSE]
    bigh <- bigh[ , 1:ncol(bigh) > nfix, drop = FALSE]
    bigh <- bigh + diag(ncol(bigh))
    baz <- diag(chol(bigh))
    return(tout$value + sum(log(baz)))
}

### like above except now t(Z) W Z is fixed matrix supplied in argument zwz
### should be symmetric and positive semidefinite
### suitable zwz can be obtained from makezwz

makezwz <- function(sigma, parm, fixed, random, obj, y, origin)
{
    checkargs(sigma, parm, fixed, random, obj, y, origin)
    nfix <- ncol(fixed)
    if (! is.list(random)) random <- list(random)
    nrand <- sapply(random, ncol)
    modmat <- fixed
    for (i in seq(along = random))
        modmat <- cbind(modmat, random[[i]], deparse.level = 0)
    scalevec <- rep(c(1, sigma), times = c(nfix, nrand))
    if (missing(y)) y <- obj$x
    bar <- mlogl(scalevec * parm, obj$pred, obj$fam, y,
        obj$root, modmat, deriv = 2, famlist = obj$famlist, origin = origin)
    zwz <- bar$hessian
    zwz <- zwz[1:nrow(zwz) > nfix, , drop = FALSE]
    zwz <- zwz[ , 1:ncol(zwz) > nfix, drop = FALSE]
    return(zwz)
}

pickle1 <- function(sigma, parm, fixed, random, obj, y, origin,
    cache, zwz, deriv = 0, ...)
{
    checkargs3part1(fixed, random, obj, y, origin)
    nfix <- ncol(fixed)
    if (! is.list(random)) random <- list(random)
    nrand <- sapply(random, ncol)
    if (! missing(cache)) {
        stopifnot(is.environment(cache))
        if (exists("parm", envir = cache, inherits = FALSE))
            parm <- cache$parm
    }
    stopifnot(is.vector(parm))
    if(length(parm) != nfix + sum(nrand))
        stop("length(parm) != number of fixed effects + number of random effects")
    idx <- seq(along = parm)
    alpha <- parm[idx <= nfix]
    cee <- parm[idx > nfix]
    checkargs3part2(alpha, cee, sigma, nfix, nrand)
    stopifnot(is.matrix(zwz))
    stopifnot(is.numeric(zwz))
    stopifnot(is.finite(zwz))
    if (any(dim(zwz) != sum(nrand)))
        stop("zwz not square matrix with dimension = number of random effects")
    stopifnot(length(deriv) == 1)
    stopifnot(deriv %in% c(0, 1))

    if (missing(y)) y <- obj$x
    trustargs <- list(objfun = penmlogl, parinit = parm, sigma = sigma,
        fixed = fixed, random = random, obj = obj, y = y, ...)
    if (is.null(trustargs$rinit)) trustargs$rinit <- 1
    if (is.null(trustargs$rmax)) trustargs$rmax <- 10
    if (! missing(origin)) trustargs$origin <- origin
    trustargs$cache <- NULL
    tout <- do.call(trust, trustargs)
    stopifnot(tout$converged)
    if ((! missing(cache))) cache$parm <- tout$argument

    alpha <- tout$argument[idx <= nfix]
    cee <- tout$argument[idx > nfix]
    pout <- pickleHelper(alpha, cee, sigma, fixed, random, obj, y,
        origin, zwz, deriv)
    if (deriv == 0)
        return(pout)
    return(list(value = pout$value, gradient = pout$gradient$pt))
}

# note: we have implementation of deriv = 2, but it doesn't work because
# of inexactness of computer arithmetic.
# thus we only allow deriv = 0 or deriv = 1

pickle2 <- function(alphasigma, parm, fixed, random, obj, y, origin,
    cache, zwz, deriv = 0, ...)
{
    checkargs3part1(fixed, random, obj, y, origin)
    nfix <- ncol(fixed)
    if (! is.list(random)) random <- list(random)
    nrand <- sapply(random, ncol)
    if (! missing(cache)) {
        stopifnot(is.environment(cache))
        if (exists("parm", envir = cache, inherits = FALSE))
            parm <- cache$parm
    }
    stopifnot(is.vector(alphasigma))
    if(length(alphasigma) != nfix + length(nrand))
        stop("length(alphasigma) != number of fixed effects + number of variance components")
    idx <- seq(along = alphasigma)
    alpha <- alphasigma[idx <= nfix]
    cee <- parm
    sigma <- alphasigma[idx > nfix]
    checkargs3part2(alpha, cee, sigma, nfix, nrand)
    stopifnot(is.matrix(zwz))
    stopifnot(is.numeric(zwz))
    stopifnot(is.finite(zwz))
    if (any(dim(zwz) != sum(nrand)))
        stop("zwz not square matrix with dimension = number of random effects")
    stopifnot(length(deriv) == 1)
    stopifnot(deriv %in% c(0, 1))

    if (missing(y)) y <- obj$x
    trustargs <- list(objfun = penmlogl2, parinit = parm, alpha = alpha,
        sigma = sigma, fixed = fixed, random = random, obj = obj, y = y, ...)
    if (is.null(trustargs$rinit)) trustargs$rinit <- 1
    if (is.null(trustargs$rmax)) trustargs$rmax <- 10
    if (! missing(origin)) trustargs$origin <- origin
    trustargs$cache <- NULL
    tout <- do.call(trust, trustargs)
    stopifnot(tout$converged)
    if ((! missing(cache))) cache$parm <- tout$argument

    cee <- tout$argument
    pout <- pickleHelper(alpha, cee, sigma, fixed, random, obj, y,
        origin, zwz, deriv)
    if (deriv == 0)
        return(pout)
    grad <- c(pout$gradient$pa, pout$gradient$pt)
    if (deriv == 1)
        return(list(value = pout$value, gradient = grad))
    paa <- pout$hessian$paa
    pac <- pout$hessian$pac
    pat <- pout$hessian$pat
    pcc <- pout$hessian$pcc
    pct <- pout$hessian$pct
    ptt <- pout$hessian$ptt
    pcc.inv <- chol2inv(chol(pcc))
    qaa <- paa - pac %*% pcc.inv %*% t(pac)
    qat <- pat - pac %*% pcc.inv %*% pct
    qtt <- ptt - t(pct) %*% pcc.inv %*% pct
    hess <- rbind(cbind(qaa, qat), cbind(t(qat), qtt))
    return(list(value = pout$value, gradient = grad, hessian = hess))
}

