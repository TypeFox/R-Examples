# PF, 2008-09-18
# is used in pfa1.R which
# computes principal factor analysis for compositional data
# Uniquenesses are nor longer of diagonal form


factanal.fit.principal1 <-
function (cmat, factors, p = ncol(cmat), start = NULL, iter.max = 10, 
    unique.tol = 1e-04) 
{
    dof <- 0.5 * ((p - factors)^2 - p - factors)
    if (dof < 0) 
        warning("negative degrees of freedom")
    if (any(abs(diag(cmat) - 1) > .Machine$single.eps)) 
        stop("must have correlation matrix")
    if (length(start)) {
        if (length(start) != p) 
            stop("start is the wrong length")
        if (any(start < 0 | start >= 1)) 
            stop("all values in start must be between 0 and 1")
        oldcomm <- 1 - start
    }
    else {
        diag(cmat) <- NA
        oldcomm <- apply(abs(cmat), 1, max, na.rm = TRUE)
    }

    # PF 10.9.2008
    H <- diag(p)-matrix(1,p,p)/p
    psi <- 1-oldcomm
    psistar <- H%*%diag(psi)%*%H
    cmatstar <- cmat-psistar
    if (iter.max < 0) 
        stop("bad value for iter.max")
    ones <- rep(1, factors)
    if (iter.max == 0) {
        z <- eigen(cmatstar, symmetric = TRUE)
        kvals <- z$values[1:factors]
        if (any(kvals <= 0)) 
            stop("impermissible estimate reached")
        Lambda <- z$vectors[, 1:factors, drop = FALSE] * rep(kvals^0.5, 
            rep.int(p, factors))
        psinew <- diag(cmat) - Lambda^2 %*% ones
    }
    if (iter.max > 0) {
        for (i in 1:iter.max) {
            z <- eigen(cmatstar, symmetric = TRUE)
            kvals <- z$values[1:factors]
            if (any(kvals <= 0)) 
                stop("impermissible estimate reached")
            Lambda <- z$vectors[, 1:factors, drop = FALSE] * 
                rep(kvals^0.5, rep.int(p, factors))
            psinew <- drop(diag(cmat) - Lambda^2 %*% ones)
            psinewstar <- H%*%diag(psinew)%*%H
            if (all(abs(psinew - psi) < unique.tol)) {
                iter.max <- i
                break
            }
            psistar <- psinewstar
            cmatstar <- cmat-psistar
        }
    }
    dn <- dimnames(cmat)[[1]]
    dimnames(Lambda) <- list(dn, paste("Factor", 1:factors, sep = ""))
    diag(cmat) <- 1
    uniq <- diag(psistar)
    names(uniq) <- dn
    ans <- list(loadings = Lambda, uniquenesses = uniq, correlation = cmat, 
        criteria = c(iterations = iter.max), factors = factors, 
        dof = dof, method = "principal",psi = psistar)
    class(ans) <- "factanal"
    ans
}
