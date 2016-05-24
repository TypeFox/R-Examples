
checkargs3part1 <- function(fixed, random, obj, y, origin)
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
    stopifnot(is.matrix(fixed))
    stopifnot(is.numeric(fixed))
    stopifnot(is.finite(fixed))
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
}

checkargs3part2 <- function(alpha, cee, sigma, nfix, nrand)
{
    if (! is.numeric(alpha))
        stop("vector of fixed effects not numeric")
    if (! all(is.finite(alpha)))
        stop("vector of fixed effects not all finite")
    if (length(alpha) != nfix)
        stop("vector of fixed effects wrong length")
    if (! is.numeric(cee))
        stop("vector of rescaled random effects not numeric")
    if (! all(is.finite(cee)))
        stop("vector of rescaled random effects not all finite")
    if (length(cee) != sum(nrand))
        stop("vector of rescaled random effects wrong length")
    if (! is.numeric(sigma))
        stop("vector of square roots of variance components not numeric")
    if (! all(is.finite(sigma)))
        stop("vector of square roots of variance components not all finite")
    if (length(sigma) != length(nrand))
        stop("vector of square roots of variance components wrong length")
}

pickleHelper <- function(alpha, cee, sigma, fixed, random, obj, y,
    origin, zwz, deriv = 0)
{
    nfix <- ncol(fixed)
    if (! is.list(random)) random <- list(random)
    nrand <- sapply(random, ncol)

    if (missing(y)) y <- obj$x

    a <- as.vector(rep(sigma, times = nrand))
    bee <- a * cee

    modmat <- fixed
    for (i in seq(along = random))
        modmat <- cbind(modmat, random[[i]], deparse.level = 0)
    mout <- mlogl(c(alpha, bee), obj$pred, obj$fam, y, obj$root,
        modmat, deriv = 2, famlist = obj$famlist, origin = origin)

    bigh <- zwz * outer(a, a) + diag(length(a))
    bigh.chol <- chol(bigh)
    val <- mout$value + sum(cee^2) / 2 + sum(log(diag(bigh.chol)))
    if (deriv == 0) return(list(value = val))

    is.alpha <- seq(along = mout$gradient) <= nfix
    is.c <- seq(along = mout$gradient) > nfix

    pa <- mout$gradient[is.alpha]
    ### Z^T (y - mu^*)
    zymoo <- mout$gradient[is.c]
    pc <- zymoo * a + cee

    bigh.inv <- chol2inv(bigh.chol)
    idx <- rep(seq(along = sigma), times = nrand)
    pt <- rep(NaN, length(sigma))
    for (k in seq(along = sigma)) {
        eek <- as.numeric(idx == k)
        pt[k] <- sum(bigh.inv * zwz * outer(a, eek)) + sum(zymoo * eek * cee)
    }
    grad <- list(pa = pa, pc = pc, pt = pt)
    if (deriv == 1) return(list(value = val, gradient = grad))

    ### M^T W^* M and M^T W^* Z
    foo <- mout$hessian[is.alpha, , drop = FALSE]
    paa <- foo[ , is.alpha, drop = FALSE]
    mwz <- foo[ , is.c, drop = FALSE]
    ### M^T W^* Z A
    pac <- sweep(mwz, 2, a, "*")
    ### Z^T W^* Z
    foo <- mout$hessian[is.c, , drop = FALSE]
    zwz.star <- foo[ , is.c, drop = FALSE]

    pcc <- zwz.star * outer(a, a) + diag(length(a))
    pat <- matrix(NaN, length(alpha), length(sigma))
    pct <- matrix(NaN, length(cee), length(sigma))
    ptt <- matrix(NaN, length(sigma), length(sigma))
    for (k in seq(along = sigma)) {
        eek <- as.numeric(idx == k)
        pat[ , k] <- as.numeric(mwz %*% cbind(eek * cee))
        pct[ , k] <- as.numeric(zwz.star %*% cbind(eek * cee)) * a - zymoo * eek
        for (m in seq(along = sigma))
            if (k <= m) {
                eem <- as.numeric(idx == m)
                qux <- sum(zwz.star * outer(eek * cee, eem * cee))
                quux <- sum(bigh.inv * zwz * outer(eek, eem))
                fook <- zwz * outer(eek, a)
                foom <- zwz * outer(eem, a)
                fook <- fook + t(fook)
                quuux <- bigh.inv %*% fook %*% bigh.inv
                quuux <- sum(quuux * foom)
                ptt[k, m] <- qux + quux - quuux
                ptt[m, k] <- qux + quux - quuux
            }
    }
    hess <- list(paa = paa, pac = pac, pat = pat, pcc = pcc,
        pct = pct, ptt = ptt)
    return(list(value = val, gradient = grad, hessian = hess))
}

pickle3 <- function(alphaceesigma, fixed, random, obj, y, origin,
    zwz, deriv = 0)
{
    checkargs3part1(fixed, random, obj, y, origin)
    nfix <- ncol(fixed)
    if (! is.list(random)) random <- list(random)
    nrand <- sapply(random, ncol)
    stopifnot(is.vector(alphaceesigma))
    if(length(alphaceesigma) != nfix + sum(nrand) + length(nrand))
        stop("length(alphaceesigma) != number of fixed effects + number of random effects + number of variance components")
    idx <- seq(along = alphaceesigma)
    alpha <- alphaceesigma[idx <= nfix]
    cee <- alphaceesigma[idx > nfix & idx <= nfix + sum(nrand)]
    sigma <- alphaceesigma[idx > nfix + sum(nrand)]
    checkargs3part2(alpha, cee, sigma, nfix, nrand)
    stopifnot(is.matrix(zwz))
    stopifnot(is.numeric(zwz))
    stopifnot(is.finite(zwz))
    if (any(dim(zwz) != sum(nrand)))
        stop("zwz not square matrix with dimension = number of random effects")
    stopifnot(length(deriv) == 1)
    stopifnot(deriv %in% c(0, 1, 2))

    pout <- pickleHelper(alpha, cee, sigma, fixed, random, obj, y,
        origin, zwz, deriv)
    if (deriv == 0)
        return(pout)
    grad <- c(pout$gradient$pa, pout$gradient$pc, pout$gradient$pt)
    if (deriv == 1)
        return(list(value = pout$value, gradient = grad))
    hess <- rbind(cbind(pout$hessian$paa, pout$hessian$pac, pout$hessian$pat),
        cbind(t(pout$hessian$pac), pout$hessian$pcc, pout$hessian$pct),
        cbind(t(pout$hessian$pat), t(pout$hessian$pct), pout$hessian$ptt))
    return(list(value = pout$value, gradient = grad, hessian = hess))
}

