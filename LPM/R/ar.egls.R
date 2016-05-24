ar.egls <-
function (x, R, order.max, na.action = na.fail, series = NULL, 
    ...) 
{
    R = sync(R)
    aic = F
    demean = T
    intercept = F
    n <- nrow(x)
    k <- ncol(x)
    if (is.null(series)) 
        series <- deparse(substitute(x))
    rescale <- TRUE
    ists <- is.ts(x)
    x <- na.action(as.ts(x))
    xfreq <- frequency(x)
    if (any(is.na(x))) 
        stop("NAs in x")
    if (ists) 
        xtsp <- tsp(x)
    x <- as.matrix(x)
    if (!is.numeric(x)) 
        stop("`x' must be numeric")
    n.used <- nrow(x)
    nser <- ncol(x)
    if (rescale) {
        sc <- sqrt(drop(apply(x, 2, var)))
        x <- x/rep(sc, rep(n.used, nser))
    }
    else sc <- rep(1, nser)
    order.max <- if (is.null(order.max)) 
        floor(10 * log10(n.used))
    else round(order.max)
    if (order.max < 0) 
        stop("order.max must be >= 0")
    if (aic) 
        order.min <- 0
    else order.min <- order.max
    A <- vector("list", order.max - order.min + 1)
    varE <- vector("list", order.max - order.min + 1)
    seA <- vector("list", order.max - order.min + 1)
    aic <- rep(Inf, order.max - order.min + 1)
    det <- function(x) {
        prod(diag(qr(x)$qr)) * (-1)^(ncol(x) - 1)
    }
    if (demean) {
        xm <- colMeans(x)
        x <- sweep(x, 2, xm)
    }
    else xm <- rep(0, nser)
    for (m in order.min:order.max) {
        y <- embed(x, m + 1)
        if (intercept) {
            if (m > 0) 
                X <- cbind(rep(1, nrow(y)), y[, (nser + 1):ncol(y)])
            else X <- as.matrix(rep(1, nrow(y)))
        }
        else {
            if (m > 0) 
                X <- y[, (nser + 1):ncol(y)]
            else X <- matrix(0, nrow(y), 0)
        }
        Y <- t(y[, 1:nser])
        N <- ncol(Y)
        XX <- t(X) %*% X
        rank <- qr(XX)$rank
        if (rank != nrow(XX)) {
            warning(paste("Model order", m))
            warning("Singularities in the computation of the projection matrix")
            warning(paste("Results are only valid up to model order", 
                m - 1))
            break
        }
        P <- if (ncol(XX) > 0) 
            solve(XX)
        else XX
        A[[m - order.min + 1]] <- Y %*% X %*% P
        YH <- A[[m - order.min + 1]] %*% t(X)
        E <- (Y - YH)
        varE[[m - order.min + 1]] <- E %*% t(E)/N
        aic[m - order.min + 1] <- n.used * log(det(varE[[m - 
            order.min + 1]])) + 2 * nser * (nser * m + intercept)
    }
    m <- which(aic == min(aic)) + order.min - 1
    y <- embed(x, m + 1)
    sigmau <- (1/(n - k * order.max - 1)) * E %*% t(E)
    for (i in 1:2) {
        T1 <- ginv(t(R) %*% (kronecker((XX), ginv(sigmau)) %*% 
            R))
        T2 <- t(R) %*% (kronecker(id(k * order.max), ginv(sigmau)))
        gam <- vec(A[[m - order.min + 1]]) + T1 %*% T2 %*% vec(E %*% 
            (X))
        beta <- R %*% gam
        AA <- unvec(beta, k, order.max)
        YH <- AA %*% t(X)
        E <- (Y - YH)
    }
    sig <- E %*% t(E)/N
    varA <- R %*% ginv(t(R) %*% (kronecker((XX/n), ginv(sig))) %*% 
        R) %*% t(R)/n
    seA[[m - order.min + 1]] <- if (ncol(varA) > 0) 
        sqrt(diag(varA))
    else numeric(0)
    if (intercept) {
        xint <- AA[, 1]
        ar <- AA[, -1]
        if (m > 0) 
            X <- cbind(rep(1, nrow(y)), y[, (nser + 1):ncol(y)])
        else X <- as.matrix(rep(1, nrow(y)))
    }
    else {
        if (m > 0) 
            X <- y[, (nser + 1):ncol(y)]
        else X <- matrix(0, nrow(y), 0)
        xint <- NULL
        ar <- AA
    }
    Y <- t(y[, 1:nser, drop = FALSE])
    YH <- AA %*% t(X)
    E <- drop(rbind(matrix(NA, m, nser), t(Y - YH)))
    aic <- aic - min(aic)
    names(aic) <- order.min:order.max
    dim(ar) <- c(nser, nser, m)
    ar <- aperm(ar, c(3, 1, 2))
    ses <- seA[[m - order.min + 1]]
    if (intercept) {
        sem <- ses[1:nser]
        ses <- ses[-(1:nser)]
    }
    else sem <- rep(0, nser)
    dim(ses) <- c(nser, nser, m)
    ses <- aperm(ses, c(3, 1, 2))
    var.pred <- varE[[m - order.min + 1]]
    if (nser > 1) {
        snames <- colnames(x)
        dimnames(ses) <- dimnames(ar) <- list(seq(length = m), 
            snames, snames)
        dimnames(var.pred) <- list(snames, snames)
        names(sem) <- colnames(E)
        snames <- colnames(E)
    }
    if (ists) {
        attr(E, "tsp") <- xtsp
        attr(E, "class") <- "ts"
    }
    if (rescale) {
        xm <- xm * sc
        if (!is.null(xint)) 
            xint <- xint * sc
        aa <- outer(sc, 1/sc)
        if (nser > 1 && m > 0) 
            for (i in 1:m) ar[i, , ] <- ar[i, , ] * aa
        var.pred <- var.pred * outer(sc, sc)
        E <- E * rep(sc, rep(NROW(E), nser))
        sem <- sem * sc
        if (m > 0) 
            for (i in 1:m) ses[i, , ] <- ses[i, , ] * aa
    }
    res <- list(order = m, ar = ar, var.pred = var.pred, x.mean = xm, 
        x.intercept = xint, aic = aic, n.used = n.used, order.max = order.max, 
        partialacf = NULL, resid = E, method = "EGLS", series = series, 
        frequency = xfreq, call = match.call(), asy.se.coef = list(x.mean = sem, 
            ar = drop(ses)))
    class(res) <- "ar"
    res
}
