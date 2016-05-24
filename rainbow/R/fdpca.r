fdpca = function (x, y, order = 2, ngrid = 500, method = "rapca", mean = TRUE, 
    level = FALSE, lambda = 3, iter = 1, ...) 
{
    n <- ncol(y)
    m <- length(x)
    if (lambda < 1) 
        stop("Lambda too small")
    if (iter < 1) 
        stop("Need at least one iteration")
    if (order < 0) 
        stop("Order must be at least 0")
    if (ngrid < n) 
        stop("Grid should be larger than number of observations per time period.")
    if (m != nrow(y)) 
        stop("x and y of incompatible dimension")
    if (order > n/2 & method != "classical") {
        warning("Not enough data for robust PCA")
        method <- "classical"
    }
    yy <- matrix(NA, nrow = ngrid, ncol = n)
    xx <- seq(min(x), max(x), l = ngrid)
    delta <- xx[2] - xx[1]
    for (i in 1:n) {
        miss <- is.na(y[, i])
        yy[, i] <- spline(x[!miss], y[!miss, i], n = ngrid)$y
    }
    if (mean) {
        if (method == "M" | method == "rapca") 
            ax <- L1median2(t(yy), method = "hoss")
        else ax <- rowMeans(yy, na.rm = TRUE)
        yy <- sweep(yy, 1, ax)
        axse <- approx(xx, sqrt(apply(yy, 1, var)/n), xout = x)$y
    }
    else axse <- NULL
    if (level) {
        bx <- colMeans(yy, na.rm = TRUE)
        yy <- sweep(yy, 2, bx)
    }
    if (level) {
        coeff <- as.matrix(bx)
        basis <- as.matrix(rep(1, m))
        colnames(coeff) <- colnames(basis) <- "level"
    }
    else coeff <- basis <- NULL
    if (mean) {
        coeff <- cbind(rep(1, n), coeff)
        basis <- cbind(approx(xx, ax, xout = x)$y, basis)
        colnames(coeff)[1] <- colnames(basis)[1] <- "mean"
    }
    if (order == 0) 
        return(list(basis = basis, coeff = coeff, weights = rep(1, 
            n), v = rep(1, n), mean.se = axse))
    if (method == "M" | method == "rapca") {
        robusteig <- rapca(yy, order = order, mean = FALSE, ...)
        B <- robusteig$coeff
        Phi <- robusteig$basis
        yyhat <- Phi %*% t(B)
        v <- colSums((yy - yyhat)^2) * delta
    }
    else {
        v <- rep(1, n)
        iter <- 1
    }
    if (method == "rapca") {
        varprop <- rep(NA, order)
        w <- rep(1, order)
    }
    else {
        for (i in 1:iter) {
            medv <- median(v)
            w <- ifelse(v < medv + lambda * sqrt(medv), 1, 0)
            V <- repmat(w/mean(w), 1, ngrid) * t(yy)
            s <- La.svd(V)
            Phi <- as.matrix(t(s$vt)[, 1:order])
        }
        varprop <- s$d^2
        varprop <- varprop/sum(s$d^2)
    }
    Phinorm = matrix(NA, length(x), order)
    Phinormngrid = matrix(NA, ngrid, order)
    for (i in 1:order) {
        Phinorm[, i] = approx(xx, Phi[, i], xout = x)$y/delta/(sqrt(sum((approx(xx, 
            Phi[, i], xout = x)$y/delta)^2)))
        Phinormngrid[, i] = approx(x, Phinorm[, i], xout = xx)$y
    }
    B <- t(yy) %*% Phinormngrid
    v <- colSums((yy - Phinormngrid %*% t(B))^2) * delta
    colnames(B) <- paste("beta", 1:order, sep = "")
    coeffdummy = B * delta
    colmeanrm = matrix(colMeans(coeffdummy), dim(B)[2], 1)
    coeff <- cbind(coeff, sweep(coeffdummy, 2, colmeanrm))
    m <- ncol(basis)
    basis = basis + Phinorm %*% colmeanrm
    for (i in 1:order) {
        basis <- cbind(basis, Phinorm[, i])
        if (sum(basis[, i + m]) < 0) {
            basis[, i + m] <- -basis[, i + m]
            coeff[, i + m] <- -coeff[, i + m]
        }
    }
    colnames(basis)[m + (1:order)] <- paste("phi", 1:order, sep = "")
    varprop <- varprop[1:order]
    yy <- yy - Phinormngrid %*% t(B)
    s <- try(La.svd(t(yy)), silent = TRUE)
    if (class(s) == "try-error") {
        s <- svd(t(yy), LINPACK = TRUE)
        s$vt <- t(s$v)
    }
    Phi2 <- as.matrix(t(s$vt)[, s$d > 1e-06])
    m <- ncol(Phi2)
    basis2 <- coeff2 <- NULL
    Phinorm2 = matrix(NA, length(x), m)
    Phinorm2ngrid = matrix(NA, ngrid, m)
    if (m > 0) {
        for (i in 1:m) {
            Phinorm2[, i] = approx(xx, Phi2[, i], xout = x)$y/delta/(sqrt(sum((approx(xx, 
                Phi2[, i], xout = x)$y/delta)^2)))
            Phinorm2ngrid[, i] = approx(x, Phinorm2[, i], xout = xx)$y
        }
        B2 <- t(yy) %*% Phinorm2ngrid
        colnames(B2) <- paste("beta", order + (1:ncol(B2)), sep = "")
        coeff2dummy <- B2 * delta
        colmeanrm2 = matrix(colMeans(coeff2dummy), dim(B2)[2], 
            1)
        coeff2 = sweep(coeff2dummy, 2, colmeanrm2)
        for (i in 1:m) {
            basis2 <- cbind(basis2, Phinorm2[, i])
            if (sum(basis2[, i]) < 0) {
                basis2[, i] <- -basis2[, i]
                coeff2[, i] <- -coeff2[, i]
            }
        }
        colnames(basis2) <- paste("phi", order + (1:m), sep = "")
    }
    return(list(basis = basis, coeff = coeff, varprop = varprop, 
        weights = w/mean(w), v = v, basis2 = basis2, coeff2 = coeff2, 
        mean.se = axse))
}

