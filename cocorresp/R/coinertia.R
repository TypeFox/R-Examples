`coinertia` <- function(y, ...)
    UseMethod("coinertia")

`coinertia.default` <- function(y, x, n.axes = NULL, weights = NULL,
                                symmetric = FALSE, ...) {
    ## weights is R0
    rsumy <- rowSums(y)
    csumy <- colSums(y)
    toty <- sum(rsumy)
    rsumx <- rowSums(x)
    csumx <- colSums(x)
    totx <- sum(rsumx)
    ## some sanity checks
    if(any(rsumy <= 0 ))
        stop("all row sums must be >0 in data matrix y")
    if(any(csumy <= 0 )) {
        y <- y[, csumy > 0, drop = FALSE]
        message("some species contain no data and were removed from data matrix y\n")
        csumy <- csumy[csumy > 0]       # colSums(y)
    }
    if(any(rsumx <= 0 ))
        stop("all row sums must be >0 in data matrix x")
    if(any(csumx <= 0 )) {
        x <- x[, csumx > 0, drop = FALSE]
        message("some species contain no data and were removed from data matrix x\n")
        csumx <- csumx[csumx > 0]       # colSums(x)
    }
    sitesy <- rownames(y)
    sitesx <- rownames(x)
    sppy <- colnames(y)
    sppx <- colnames(x)
    y <- data.matrix(y)
    x <- data.matrix(x)
    nrx <- nrow(y)
    nry <- nrow(x)
    ncy <- ncol(y)
    ncx <- ncol(x)
    if (nrx != nry) {
        stop("Number of rows in y and x is not equal")
    }
    max.axes <- min(ncy, ncx, nry, nrx) - 1
    if (is.null(n.axes)) {
        n.axes <- max.axes
    } else {
        if (n.axes > max.axes) {
            n.axes <- max.axes
            warning("n.axes greater than min(n,p,q)-1,\nreset to min(n,p,q)-1")
        }
    }
    Axes <- seq_len(n.axes)
    ax.names <- paste("COIN", Axes, sep = "")
    if (is.null(weights)) {
        if (symmetric) {
            weights <- (rsumy + rsumx) / 2
        } else {
            weights <- rsumy
        }
    }
    .R0 <- weights / sum(weights)
    .csy <- csumy / toty
    .csx <- csumx / totx
    Q1 <- diag(toty / rsumy) %*% y %*% diag(1 / csumy) - 1
    Q2 <- diag(totx / rsumx) %*% x %*% diag(1 / csumx) - 1
    colnames(Q1) <- sppy
    colnames(Q2) <- sppx
    rownames(Q1) <- sitesx
    rownames(Q2) <- sitesy
    rooty <- sqrt(.csy)
    rootx <- sqrt(.csx)
    Dy <- diag(rooty)
    Dx <- diag(rootx)
    A <- Dy %*% t(Q1) %*% diag(.R0) %*% Q2 %*% Dx
    svdA <- La.svd(A)
    U <- diag(1 / rooty) %*% svdA$u
    V <- diag(1 / rootx) %*% t(svdA$vt)
    Ksi <- Q1 %*% diag(.csy) %*% U
    Psi <- Q2 %*% diag(.csx) %*% V
    L <- svdA$d^2
    U1 <- U[, Axes, drop = FALSE]
    U2 <- V[, Axes, drop = FALSE]
    colnames(U1) <- colnames(U2) <- ax.names
    rownames(U1) <- colnames(Q1)
    rownames(U2) <- colnames(Q2)
    X1 <- Ksi[, Axes, drop = FALSE]
    X2 <- Psi[, Axes, drop = FALSE]
    colnames(X1) <- colnames(X2) <- ax.names
    lambda <- L[Axes]
    names(lambda) <- ax.names
    res <- list(scores = list(species = list(Y = U1, X = U2),
                sites = list(Y = X1, X = X2)), weights = weights,
                lambda = lambda, n.axes = n.axes,
                symmetric = symmetric, call = match.call())
    class(res) <- c("coinertia", "list")
    res
}

`print.coinertia` <- function(x, digits = max(3, getOption("digits") - 3),
                              ...) {
    writeLines("\nCoinertia Analysis\n")
    writeLines(strwrap(pasteCall(x$call)))
    cat("\nEigenvalues:\n")
    print(round(eigenvals(x), digits), ..., print.gap = 3)
}

`eigenvals.coinertia` <- function(x, ...) {
    x$lambda
}
