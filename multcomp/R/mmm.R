
### value added version of dlm:bdiag by Giovanni Petris <GPetris@uark.edu>
.bdiag <- function (...)  {
    if (nargs() == 1) 
        x <- as.list(...)
    else x <- list(...)
    n <- length(x)
    if (n == 0) 
        return(NULL)
    x <- lapply(x, function(y) if (length(y)) 
        as.matrix(y)
    else stop("Zero-length component in x"))
    d <- array(unlist(lapply(x, dim)), c(2, n))
    rr <- d[1, ]
    cc <- d[2, ]
    rsum <- sum(rr)
    csum <- sum(cc)
    out <- array(0, c(rsum, csum))
    ind <- array(0, c(4, n))
    rcum <- cumsum(rr)
    ccum <- cumsum(cc)
    ind[1, -1] <- rcum[-n]
    ind[2, ] <- rcum
    ind[3, -1] <- ccum[-n]
    ind[4, ] <- ccum
    imat <- array(1:(rsum * csum), c(rsum, csum))
    iuse <- apply(ind, 2, function(y, imat) imat[(y[1] + 1):y[2], 
        (y[3] + 1):y[4]], imat = imat)
    iuse <- as.vector(unlist(iuse))
    out[iuse] <- unlist(x)
    rownames(out) <- unlist(sapply(x, rownames))
    colnames(out) <- unlist(sapply(x, colnames))
    return(out)
}

### collect multiple marginal models
mmm <- function(...) {

     ret <- list(...)
     if (is.null(names(ret))) 
         names(ret) <- as.character(match.call(expand.dots = TRUE))[-1]
     class(ret) <- "mmm"
     ret
}

### collect multiple linear functions
mlf <- function(...) {

     ret <- list(...)
     class(ret) <- "mlf"
     ret
}

### extract coefs
coef.mmm <- function(object, ...) {

     ret <- lapply(object, function(o) modelparm(o)$coef)
     n <- lapply(1:length(ret), function(i) 
                 paste(names(object)[i], names(ret[[i]]), sep = ": "))
     ret <- unlist(ret)
     names(ret) <- unlist(n)
     ret
}

Estfun <- function(x, ...) {
     e <- estfun(x, ...)
     e[is.na(e)] <- 0
     e
}

### extract estimating functions
estfun.mmm <- function(x, ...)
    do.call("cbind", lapply(x, Estfun, ...))

### extract bread
bread.mmm <- function(x, ...)
    .bdiag(lapply(x, bread))

### set-up total covariance matrix
vcov.mmm <- function(object, ...) {

     ret <- sandwich(object, bread. = bread, meat. = meat, adjust = FALSE)
     D <- do.call("c", lapply(object, function(o) diag(vcov(o))))
     D <- diag(sqrt(D))
     ret <- D %*% cov2cor(ret) %*% D
     cf <- coef(object)
     rownames(ret) <- colnames(ret) <- names(cf)
     ret
}

### construct all linear functions and call glht
### Note: in glht check model for class mmm and
### call this function
### need to export and document mmm and mlf only
glht.mlf <- function(model, linfct, ...) {

    stopifnot(inherits(model, "mmm"))
    if (length(linfct) == 1) {
        linfct <- linfct[rep(1, length(model))]
        names(linfct) <- names(model)
    }
    if (!isTRUE(all.equal(names(model), names(linfct))))
        stop("names of ", sQuote("model"), " and ", sQuote("linfct"),
             " are not identical")

    K <- lapply(names(model), function(i) 
                glht(model[[i]], linfct = linfct[[i]], ...)$linfct)
    for (k in 1:length(K)) {
        rownames(K[[k]]) <- paste(names(model)[k], 
            rownames(K[[k]]), sep = ": ")
        colnames(K[[k]]) <- paste(names(model)[k], 
            colnames(K[[k]]), sep = ": ")
    }

    glht.matrix(model, linfct = .bdiag(K), ...)
}
