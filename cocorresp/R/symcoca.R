`symcoca` <- function(y, x, n.axes = NULL, R0 = NULL, symmetric = FALSE,
                      nam.dat = NULL) {
    ## Y1 is y and Y2 is x - this matters as this is a symmetric coca
    ## but the analysis is still described in favour of the response
    ## The symmetric argument sets whether the weights (R0)
    ## are the average of the rowsums for x and y or if the rowsums of y,
    ## the response,are used for R0.
    ## the role of x or y does not matter then if symmetric = TRUE
    if(is.null(nam.dat)) {
        namY <- deparse(substitute(y))
        namX <- deparse(substitute(x))
    } else {
        namY <- nam.dat$namY
        namX <- nam.dat$namX
    }
    rsum1 <- rowSums(y)
    csum1 <- colSums(y)
    tot1 <- sum(rsum1)
    rsum2 <- rowSums(x)
    csum2 <- colSums(x)
    site.names1 <- rownames(y)
    site.names2 <- rownames(x)
    spp.names1 <- colnames(y)
    spp.names2 <- colnames(x)
    y <- as.matrix(y)
    x <- as.matrix(x)
    tot2 <- sum(rsum2)
    nrow1 <- nrow(y)
    nrow2 <- nrow(x)
    if(nrow1 != nrow2) stop("Number of rows in y and x is not equal")
    max.axes <- min((p <- ncol(y)), (q <- ncol(x)), (n.row <- nrow(y)), nrow(x)) - 1
    if(is.null(n.axes)) {
        n.axes <- max.axes
    } else {
        if(n.axes > max.axes) {
            n.axes <- max.axes
            warning("n.axes greater than min(n,p,q)-1,\nreset to min(n,p,q)-1")
        }
    }
    ax.names <- paste("COCA", 1:n.axes, sep = " ")
    if(is.null(R0)) {
        if (symmetric)
            R0 <- (rsum1 + rsum2) / 2
        else R0 <- rsum1
    }
    .R0 <- R0 / sum(R0)
    .csum1 <- csum1 / tot1
    .csum2 <- csum2 / tot2
    Q1 <- diag(tot1 / rsum1) %*% y %*% diag(1 / csum1) - 1
    Q2 <- diag(tot2 / rsum2) %*% x %*% diag(1 / csum2) - 1
    colnames(Q1) <- spp.names1
    colnames(Q2) <- spp.names2
    rownames(Q1) <- site.names1
    rownames(Q2) <- site.names2
    res.coin <- fitCoinertia(Q1, .csum1, Q2, .csum2, .R0, n.axes = n.axes)
    X <- (res.coin$scores$site$Y + res.coin$scores$site$X) / 2
    load1 <- Conj(t(Q1)) %*% diag(.R0) %*% X %*% solve(Conj(t(X)) %*% diag(.R0) %*% X)
    load2 <- Conj(t(Q2)) %*% diag(.R0) %*% X %*% solve(Conj(t(X)) %*% diag(.R0) %*% X)
    rownames(load1) <- spp.names1
    rownames(load2) <- spp.names2
    colnames(load1) <- colnames(load2) <- ax.names
    Q1r <- Q1 - X %*% Conj(t(load1))
    Q2r <- Q2 - X %*% Conj(t(load2))
    tot.inK1 <- .csum1 %*% diag(Conj(t(Q1)) %*% diag(.R0) %*% Q1)
    tot.inK2 <- .csum2 %*% diag(Conj(t(Q2)) %*% diag(.R0) %*% Q2)
    tot.inertia <- list(Y = tot.inK1, X = tot.inK2)
    res.inK1 <- .csum1 %*% diag(Conj(t(Q1r)) %*% diag(.R0) %*% Q1r)
    res.inK2 <- .csum2 %*% diag(Conj(t(Q2r)) %*% diag(.R0) %*% Q2r)
    res.inertia <- list(Y = res.inK1, X = res.inK2)
    retval <- list(scores = res.coin$scores,
                   lambda = res.coin$lambda, #coinertia = res.coin,
                   X = X,
                   loadings = list(Y = load1, X = load2),
                   residuals = list(Y = Q1r, X = Q2r),
                   inertia = list(total = tot.inertia,
                   residual = res.inertia),
                   rowsum = list(rsum1 = rsum1, rsum2 = rsum2),
                   colsum = list(csum1 = csum1, csum2 = csum2),
                   nam.dat = list(namY = namY, namX = namX),
                   n.axes = n.axes,
                   weights = .R0,
                   call = match.call())
    class(retval) <- c("symcoca", "coca", "list")
    retval
}
