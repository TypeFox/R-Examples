PGDnpmr <-
function(B, b, X, Y, lambda, s, group = NULL, accelerated = TRUE,
    eps = 1e-7, maxit = 1e5, quiet = TRUE) {

    if (is.null(group)) group = rep(1, nrow(B))
    sumY = colSums(Y)
    XtY = crossprod(X, Y)
    W = which(Y == 1)
    eta = matrix(1, nrow(X), 1) %*% t(b) + X %*% B
    P = exp(eta)/rowSums(exp(eta))

    objectivePath = c(objectiveFast(B, P, W, lambda), rep(NA, maxit))
    it = 0
    diff = eps + 1
    Sys.time = Sys.time()
    while(abs(diff) > eps & it < maxit) {
        it = it + 1
        if (!quiet) print(objectivePath[it])
        B. = prox(B + s*XtY - s*crossprod(X, P), s*lambda, group)
        b. = b + s*sumY - s*colSums(P)
        if (accelerated) {
            B. = B. + it/(it+3)*(B. - B)
            b. = b. + it/(it+3)*(b. - b)
        }
        eta = t(t(as.matrix(X%*%B.)) + b.)
        expeta = exp(eta)
        P. = expeta/rowSums(expeta)
        objectivePath[it+1] = objectiveFast(B., P., W, lambda)
        while(objectivePath[it+1] > objectivePath[it]){
            if (!quiet) print(paste('s =', s))
            s = s/2
            B. = prox(B + s*XtY - s*crossprod(X, P), s*lambda, group)
            b. = b + s*sumY - s*colSums(P)
            if (accelerated) {
                B. = B. + it/(it+3)*(B. - B)
                b. = b. + it/(it+3)*(b. - b)
            }
            eta = t(t(as.matrix(X%*%B.)) + b.)
            expeta = exp(eta)
            P. = expeta/rowSums(expeta)
            objectivePath[it+1] = objectiveFast(B., P., W, lambda)
        }
        B = B.
        b = b. - mean(b.)
        P = P.
        diff = (objectivePath[it] - objectivePath[it+1])/objectivePath[it+1]
    }

    list(B = B, b = b, objectivePath = objectivePath[1:(it+1)],
        time = Sys.time() - Sys.time)
}
