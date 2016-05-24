### implements (7) of the design doc

quickle <- function(alphanu, bee, fixed, random, obj, y, origin, zwz,
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
    stopifnot(length(deriv) == 1)
    stopifnot(deriv %in% c(0, 1, 2))
    stopifnot(is.vector(alphanu))
    stopifnot(is.numeric(alphanu))
    stopifnot(is.finite(alphanu))
    if (length(alphanu) != nfix + length(nrand))
        stop("alphanu wrong length")
    stopifnot(is.vector(bee))
    stopifnot(is.numeric(bee))
    stopifnot(is.finite(bee))
    if (length(bee) != sum(nrand))
        stop("bee wrong length")

    idx <- seq(along = alphanu)
    is.alpha <- idx <= nfix
    is.nu <- nfix < idx
    alpha <- alphanu[is.alpha]
    nu <- alphanu[is.nu]
    dee <- rep(nu, times = nrand)

    if (any(nu < 0))
        return(list(value = Inf, alpha = alpha, bee = bee, nu = nu))

    modmat <- cbind(fixed, Reduce(cbind, random))
    ### note: despite documentation of the mlogl function, it actually
    ### works to have modmat a matrix rather than a 3-way array
    mout <- mlogl(c(alpha, bee), obj$pred, obj$fam, y, obj$root, modmat,
        deriv = 2, famlist = obj$famlist, origin = origin)

    idx <- seq(along = mout$gradient)
    is.alpha <- idx <= nfix
    is.bee <- nfix < idx


    if (missing(origin))
        mymlogl <- function(bee)
            mlogl(c(alpha, bee), obj$pred, obj$fam, y, obj$root,
                modmat, deriv = 2, famlist = obj$famlist)
    else
        mymlogl <- function(bee)
            mlogl(c(alpha, bee), obj$pred, obj$fam, y, obj$root,
                modmat, deriv = 2, famlist = obj$famlist, origin = origin)
    objfun <- function(bee) {
        ### note: despite documentation of the mlogl function, it actually
        ### works to have modmat a matrix rather than a 3-way array
        mout <- mymlogl(bee)
        val <- mout$value + sum(bee^2 / dee) / 2
        grad <- mout$gradient[is.bee] + bee / dee
        hess <- mout$hessian
        hess <- hess[is.bee, , drop = FALSE]
        hess <- hess[ , is.bee, drop = FALSE]
        ### see note on help page for diag !!!
        hess <- hess + diag(1 / dee, nrow = length(dee))
        return(list(value = val, gradient = grad, hessian = hess))
    }
    tout <- trust(objfun, bee, rinit = 1, rmax = 10, iterlim = 1000)
    stopifnot(tout$converged)
    bee <- tout$argument
    mout <- mlogl(c(alpha, bee), obj$pred, obj$fam, y, obj$root, modmat,
        deriv = 2, famlist = obj$famlist, origin = origin)

    a <- sqrt(dee)
    bigh <- zwz * outer(a, a) + diag(length(a))
    bigh.chol <- chol(bigh)

    val <- mout$value + sum(bee^2 / dee) / 2 + sum(log(diag(bigh.chol)))
    if (deriv == 0)
        return(list(value = val, alpha = alpha, bee = bee, nu = nu))

    pa <- mout$gradient[is.alpha]
    pb <- mout$gradient[is.bee] + bee / dee

    bigh.inv <- chol2inv(bigh.chol)
    idx <- rep(seq(along = nu), times = nrand)
    pn <- rep(NaN, length(nu))
    for (k in seq(along = nu)) {
        eek <- as.numeric(idx == k)
        pn[k] <- sum(bigh.inv * zwz * outer(a, eek / a)) / 2 -
            sum(bee^2 / dee^2 * eek) / 2
    }

    if (deriv == 1)
        return(list(value = val, gradient = c(pa, pn),
            alpha = alpha, bee = bee, nu = nu))

    foo <- mout$hessian[is.alpha, , drop = FALSE]
    paa <- foo[ , is.alpha, drop = FALSE]
    pab <- foo[ , is.bee, drop = FALSE]
    foo <- mout$hessian[is.bee, , drop = FALSE]
    pbb <- foo[ , is.bee, drop = FALSE] + diag(1 / dee, nrow = length(dee))
    pan <- matrix(0, length(alpha), length(nu))
    pbn <- matrix(NaN, length(bee), length(nu))
    pnn <- matrix(NaN, length(nu), length(nu))
    bigh.inverse <- chol2inv(bigh.chol)
    for (k in seq(along = nu)) {
        eek <- as.numeric(idx == k)
        pbn[ , k] <- (- bee * eek / dee^2)
        fook <- zwz * outer(a, eek / a)
        fook <- bigh.inverse %*% fook
        for (m in seq(along = nu)) {
            eem <- as.numeric(idx == m)
            foom <- zwz * outer(a, eem / a)
            foom <- bigh.inverse %*% foom
            pnn[k, m] <- sum(bee^2 * eek * eem / dee^3) -
                sum(fook * t(foom)) / 2
        }
    }
    pbb.inv <- chol2inv(chol(pbb))
    poo <- rbind(cbind(paa, pan), cbind(t(pan), pnn))
    pob <- rbind(pab, t(pbn))
    hess <- poo - pob %*% pbb.inv %*% t(pob)
    return(list(value = val, gradient = c(pa, pn), hessian = hess,
            alpha = alpha, bee = bee, nu = nu, pbb.inv = pbb.inv,
            pba = t(pab), pbn = pbn))
}

