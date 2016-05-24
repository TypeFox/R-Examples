MultLogLik <-
function(x, nvec, paramch=1) {

    k <- length(x)
    p <- invbaselogit(x[-k])$probs

    if (paramch==1) {
        f <- (exp(x[k])-1)/(exp(x[k])+1)
    }
    if (paramch==2) {
        minp <- min(p)
        minf <- -minp/(1-minp)
        f <- (exp(x[k])+minf)/(exp(x[k])+1)
    }

    list("MultLogLik"=MultLogLikP(p, f, nvec))
}

MultLogLikP <-
function(p, f, nvec) {

    k <- length(p)
    if (length(nvec) != k * (k+1)/2) {
        stop("length mistmatch between p and nvec")
    }
    if (is.matrix(p)) {
        if (nrow(p) == 1) p <- c(p)
        else stop("p should be a vector")
    }
    minp <- min(p)
    minf <- -minp/(1-minp)

    if (is.na(minp) || is.na(f) || f < minf) {
        MultLogLik <- -Inf
    }
    else {

        Q <- (1-f) * p %*% t(p) + f * diag(p)
        R <- rep(1:k, times=k)
        C <- rep(1:k, each=k)
        Q[R > C] <- 2 * Q[R > C]
        qvec <- Q[R >= C]

        if (any(qvec <= 0 | qvec > 1)) {
            MultLogLik <- -Inf
        }
        else {
            MultLogLik <- sum(nvec * log(qvec))
        }
    }

    return(MultLogLik)
}

