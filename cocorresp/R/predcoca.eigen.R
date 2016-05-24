"predcoca.eigen" <- function (y, x, R0 = NULL, n.axes = NULL,
                              nam.dat = NULL) {
    ## y is the response matrix
    ## x is the predictor matrix
    if (is.null(nam.dat)) {
        namY <- deparse(substitute(y))
        namX <- deparse(substitute(x))
    } else {
        namY <- nam.dat$namY
        namX <- nam.dat$namX
    }
    nam.dat <- list(namY = namY, namX = namX)
    site.namesY <- rownames(y)
    site.namesX <- rownames(x)
    spp.namesY <- colnames(y)
    spp.namesX <- colnames(x)
    y <- as.matrix(y)
    x <- as.matrix(x)
    y <- y / sum(rowSums(y))
    x <- x / sum(rowSums(x))
    Yrsum <- rowSums(y)
    Ycsum <- colSums(y)
    Ytot <- sum(Yrsum)
    Xrsum <- rowSums(x)
    Xcsum <- colSums(x)
    Xtot <- sum(Xrsum)
    if (is.null(R0)) {
        R0 <- Yrsum
    }
    Ychi1 <- mcChi(y, R0)$Ychi
    Ychi2 <- mcChi(x, R0)$Ychi
    p <- ncol(y)
    q <- ncol(x)
    n.row <- nrow(y)
    max.axes <- min(p, q, n.row, nrow(x)) - 1
    if (is.null(n.axes)) {
        n.axes <- max.axes
    } else {
        if (n.axes > max.axes) {
            n.axes <- max.axes
            warning("n.axes greater than min(n,p,q)-1, reset to min(n,p,q)-1")
        }
    }
    nq1 <- min(n.row, q) - 1
    nq1 <- min(nq1, n.axes)
    R <- Yrsum * Xrsum / R0
    invR <- diag(1 / R)
    invYc <- diag(1 / Ycsum)
    invXc <- diag(1 / Xcsum)
    Z <- Conj(t(y)) %*% invR %*% x
    C1 <- rep(1, length(Ycsum)) %*% diag(Ycsum)
    C2 <- rep(1, length(Xcsum)) %*% diag(Xcsum)
    lambda <- numeric(n.axes)
    U1 <- matrix(0, nrow = length(Ycsum), ncol = n.axes)
    U2 <- matrix(0, nrow = length(Xcsum), ncol = n.axes)
    X1 <- matrix(0, nrow = n.row, ncol = n.axes)
    X2 <- X1
    if (p < q) {
        for (s in seq_len(nq1)) {
            PI1 <- Conj(t(C1)) %*% solve(C1 %*% invYc %*% Conj(t(C1))) %*%
                C1 %*% invYc
            PI2 <- Conj(t(C2)) %*% solve(C2 %*% invXc %*% Conj(t(C2))) %*%
                C2 %*% invXc
            ZKZt <- (diag(1, nrow = nrow(PI1),
                          ncol = ncol(PI1)) - PI1 ) %*% Z %*% invXc %*%
                    (diag(1, nrow = nrow(PI2),
                          ncol = ncol(PI2)) - PI2) %*% Conj(t(Z))
            eig <- eigen(solve(diag(Ycsum),ZKZt ))
            eigVals <- diag(Re(eig$values))
            maxL <- max(eigVals)
            eigmax <- which.max(eigVals)
            u <- Re(eig$vectors[, eigmax, drop = FALSE])
            v <- invXc %*% (diag(1, nrow = nrow(PI2),
                                 ncol = ncol(PI2)) - PI2) %*% Conj(t(Z)) %*% u
            u <- u %*% diag(1 / sqrt(Conj(t(u)) %*% diag(Ycsum) %*% u))
            v <- v %*% diag(1 / sqrt(Conj(t(v)) %*% diag(Xcsum) %*% v))
            x1 <- diag(1 / Yrsum) %*% y %*% u
            x2 <- diag(1 / Xrsum) %*% x %*% v
            lambda[s] <- maxL
            U1[, s] <- u
            U2[, s] <- v
            X1[, s] <- x1
            X2[, s] <- x2
            c2 <- Conj(t(x2)) %*% diag(R0 / Xrsum) %*% x
            C2 <- rbind(C2, c2)
        }
    } else {
       for (s in seq_along(nq1)) {
           PI1 <- Conj(t(C1)) %*% solve(C1 %*% invYc %*% Conj(t(C1))) %*%
               C1 %*% invYc
           PI2 <- Conj(t(C2)) %*% solve(C2 %*% invXc %*% Conj(t(C2))) %*%
               C2 %*% invXc
           ZtKZ <- (diag(1, nrow = nrow(PI2),
                         ncol = ncol(PI2)) - PI2) %*% Conj(t(Z)) %*% invYc %*%
                   (diag(1, nrow = nrow(PI1),
                         ncol = ncol(PI1)) - PI1) %*% Z
           eig <- eigen(solve(diag(Xcsum), ZtKZ))
           eigVals <- diag(Re(eig$values))
           maxL <- max(eigVals)
           eigmax <- which.max(eigVals)
           v <- Re(eig$vectors[, eigmax, drop = FALSE])
           u <- invYc %*% (diag(1, nrow = nrow(PI1),
                                ncol = ncol(PI1)) - PI1) %*% Z %*% v
           u <- u %*% diag(1 / sqrt(Conj(t(u)) %*% diag(Ycsum) %*% u))
           v <- v %*% diag(1 / sqrt(Conj(t(v)) %*% diag(Xcsum) %*% v))
           x1 <- diag(1 / Yrsum) %*% y %*% u
           x2 <- diag(1 / Xrsum) %*% x %*% v
           lambda[s] <- maxL
           U1[, s] <- u
           U2[, s] <- v
           X1[, s] <- x1
           X2[, s] <- x2
           c2 <- Conj(t(x2)) %*% diag(R0 / Xrsum) %*% x
           C2 <- rbind(C2, c2)
       }
    }
    rownames(U1) <- spp.namesY
    rownames(U2) <- spp.namesX
    rownames(X1) <- site.namesY
    rownames(X2) <- site.namesX
    ax.names <- paste("COCA", 1:n.axes, sep = " ")
    colnames(U1) <- colnames(U2) <- colnames(X1) <- colnames(X2)<- ax.names
    names(lambda) <- ax.names
    retval <- list(nam.dat = nam.dat, method = "eigen",
                   n.axes = n.axes, lambda = lambda,
                   scores = list(species = list(Y = U1, X = U2),
                                 site = list(Y = X1, X = X2)),
                   R0 = R0, Ychi = list(Ychi1 = Ychi1, Ychi2 = Ychi2),
                   call = match.call(), loadings = NULL, fitted = NULL,
                   varianceExp = NULL, totalVar = NULL)
    class(retval) <- c("predcoca", "coca", "list")
    retval
}

