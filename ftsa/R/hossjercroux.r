hossjercroux = function (X, tol = 1e-06, maxstep = 100, na.rm = TRUE)
{
    n <- nrow(X)
    p <- ncol(X)
    m = apply(X, 2, median.default, na.rm = na.rm)
    hctol <- max(1, min(abs(m), na.rm = na.rm)) * tol
    for (k in 1:maxstep) {
        mold <- m
        XX <- sweep(X, 2, m)
        dx <- norme(XX)
        if (min(abs(dx)) > tol)
            w <- 1/dx
        else {
            w <- rep(0, n)
            w[dx > tol] <- 1/dx[dx > tol]
        }
        delta <- colSums(XX * repmat(w/sum(w), 1, p), na.rm = na.rm)
        nd <- sqrt(sum(delta^2))
        maxhalf <- ifelse(nd < hctol, 0, log2(nd/hctol))
        m <- mold + delta
        nstep <- 0
        oldmobj <- mrobj(X, mold)
        while ((mrobj(X, m) > oldmobj) & (nstep <= maxhalf)) {
            nstep <- nstep + 1
            m <- mold + delta/(2^nstep)
        }
        if (nstep > maxhalf)
            return(mold)
    }
    return(mold)
}
