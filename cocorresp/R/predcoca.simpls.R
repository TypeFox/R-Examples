`predcoca.simpls` <- function(y, x, R0 = NULL, n.axes = NULL,
                              nam.dat = NULL) {
    if(is.null(nam.dat)) {
        namY <- deparse(substitute(y))
        namX <- deparse(substitute(x))
    } else {
        namY <- nam.dat$namY
        namX <- nam.dat$namX
    }
    nam.dat = list(namY = namY, namX = namX)
    site.names1 <- rownames(y)
    site.names2 <- rownames(x)
    spp.names1 <- colnames(y)
    spp.names2 <- colnames(x)
    y <- as.matrix(y)
    x <- as.matrix(x)
    dimnames(y) <- dimnames(x) <- NULL
    Yrsum <- rowSums(y)
    Ycsum <- colSums(y)
    Ytot <- sum(Yrsum)
    Xrsum <- rowSums(x)
    Xcsum <- colSums(x)
    Xtot <- sum(Xrsum)
    if (is.null(R0))
        R0 <- Yrsum
    p <- ncol(y)
    q <- ncol(x)
    n.row <- nrow(y)
    max.axes <- min(p, q, n.row, nrow(x)) - 1
    if(is.null(n.axes)) {
        n.axes <- max.axes
    } else {
        if(n.axes > max.axes) {
            n.axes <- max.axes
            warning("n.axes greater than min(n,p,q)-1, reset to min(n,p,q)-1")
        }
    }
    R0.scaled <- R0 / sum(R0)
    Ycsum.scaled <- Ycsum / Ytot
    Xcsum.scaled <- Xcsum / Xtot
    Ychi1 <- mcChi(y, R0)
    Ychi2 <- mcChi(x, R0)
    pls.res <- simpls(Ychi2$Ychi, Ychi1$Ychi, n.axes)
    U2 <- diag(1 / sqrt(Xcsum.scaled)) %*% pls.res$projection
    X2 <- diag(1 / Xrsum) %*% x %*% U2
    U1 <- diag(1 / sqrt(Ycsum.scaled)) %*% pls.res$Yloadings
    X1 <- diag(1 / Yrsum) %*% y %*% U1
    loadings1 <- U1
    loadings2 <- diag(1 / sqrt(Xcsum.scaled)) %*% pls.res$loadings
    Yhat <- Ychi2$Ychi %*% pls.res$coefficients[, , n.axes]
    Yhat1 <- diag(1 / sqrt(R0.scaled)) %*% Yhat %*% diag(1 / sqrt(Ycsum.scaled))
    Yhat1 <- Yhat1 + matrix(1, n.row, 1) %*% matrix(1, 1, p)
    Yhat1 <- diag(Yrsum) %*% Yhat1 %*% diag(Ycsum / Ytot)
    rownames(U1) <- spp.names1
    rownames(U2) <- spp.names2
    rownames(X1) <- site.names1
    rownames(X2) <- site.names2
    rownames(loadings1) <- spp.names1
    rownames(loadings2) <- spp.names2
    rownames(Yhat1) <- rownames(Yhat) <- site.names1
    colnames(Yhat1) <- colnames(Yhat) <- spp.names1
    retval <- list(nam.dat = nam.dat, call = match.call(), method = "simpls",
                   scores = list(species = list(Y = U1, X = U2),
                   site = list(Y = X1, X = X2)),
                   loadings = list(Y = loadings1, X = loadings2),
                   fitted = list(Yhat = Yhat, Yhat1 = Yhat1),
                   varianceExp = list(Xblock = pls.res$Xvar,
                   Yblock = pls.res$Yvar),
                   totalVar = list(Xblock = pls.res$Xtotvar,
                   Yblock = pls.res$Ytotvar),
                   lambda = NULL, n.axes = n.axes,
                   Ychi = list(Ychi1 = Ychi1$Ychi, Ychi2 = Ychi2$Ychi),
                   R0 = R0)
    class(retval) <- c("predcoca", "coca", "list")
    retval
}

