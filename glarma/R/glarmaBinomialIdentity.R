glarmaBinomialIdentity <- function(y, X, offset = NULL, delta, phiLags, thetaLags,
                                   method = "FS") {
    n.trial <- y[, 1] + y[, 2]
    ## Note this was changed from earlier versions for consistency with R's binomial
    ## glm conventions.
    r     <- ncol(X)
    y     <- y[, 1]
    n     <- length(y)
    p     <- length(phiLags)
    q     <- length(thetaLags)
    beta  <- delta[1:r]
    phi   <- delta[(r + 1):(r + p)]
    theta <- delta[(r + p + 1):(r + p + q)]
    mpq   <- 0
    if ((p + q) > 0) {
        mpq <- max(phiLags[p], thetaLags[q])
    }
    nmpq <- n + mpq
    s    <- r + p + q
    e    <- array(0, nmpq)
    Z    <- array(0, nmpq)
    W    <- array(0, nmpq)
    mu   <- array(0, nmpq)
    e.d  <- array(0, c(s, nmpq))
    Z.d  <- array(0, c(s, nmpq))
    W.d  <- array(0, c(s, nmpq))
    if (method == "NR") {
        e.dd <- array(0, c(s, s, nmpq))
        Z.dd <- array(0, c(s, s, nmpq))
        W.dd <- array(0, c(s, s, nmpq))
    }
    if(is.null(offset)) eta<-X %*% beta else eta<- X %*% beta + offset
    ll    <- 0
    ll.d  <- matrix(0, ncol = 1, nrow = s)
    ll.dd <- matrix(0, ncol = s, nrow = s)
    ## Start loop on time.
    for (time in 1:n) {
        tmpq <- time + mpq
        if (p > 0) {
            Z.d[(r + 1):(r + p), tmpq] <- Z[tmpq - phiLags] + e[tmpq - phiLags]
            if (method == "NR") {
                Z.dd[(r + 1):(r + p), , tmpq] <- t((Z.d + e.d)[, (tmpq - phiLags)])
                Z.dd[, (r + 1):(r + p), tmpq] <- Z.dd[, (r + 1):(r + p), tmpq] +
                                                 (Z.d + e.d)[, (tmpq - phiLags)]
            }
            for (i in 1:p) {
                Z[tmpq] <- Z[tmpq] + phi[i] * (Z + e)[tmpq - phiLags[i]]
                Z.d[, tmpq] <- Z.d[, tmpq] + phi[i] * (Z.d[, tmpq - phiLags[i]] +
                  e.d[, tmpq - phiLags[i]])
                if (method == "NR") {
                  Z.dd[, , tmpq] <- Z.dd[, , tmpq] + phi[i] *
                                    (Z.dd[, , tmpq - phiLags[i]] +
                                    e.dd[, , tmpq - phiLags[i]])
                }
            }
        }
        if (q > 0) {
            Z.d[(r + p + 1):(r + p + q), tmpq] <- e[tmpq - thetaLags]
            if (method == "NR") {
                Z.dd[(r + p + 1):(r + p + q), , tmpq] <- Z.dd[(r + p + 1):
                                                         (r + p + q), , tmpq] +
                                                         t(e.d[, tmpq - thetaLags])
                Z.dd[, (r + p + 1):(r + p + q), tmpq] <- Z.dd[, (r + p + 1):
                                                         (r + p + q), tmpq] +
                                                         e.d[, tmpq - thetaLags]
            }
            for (i in 1:q) {

                Z[tmpq] <- Z[tmpq] + theta[i] * e[tmpq - thetaLags[i]]
                Z.d[, tmpq] <- Z.d[, tmpq] + theta[i] * e.d[, tmpq - thetaLags[i]]
                if (method == "NR") {
                  Z.dd[, , tmpq] <- Z.dd[, , tmpq] + theta[i] *
                                    e.dd[, , tmpq - thetaLags[i]]
                }
            }
        }

        W[tmpq] <- eta[time] + Z[tmpq]

        W.d[, tmpq] <- matrix(c(X[time, ], rep(0, p + q)), ncol = 1) + Z.d[, tmpq]
        if (method == "NR") {

            W.dd[, , tmpq] <- Z.dd[, , tmpq]
        }
        pt <- exp(W[tmpq])/(1 + exp(W[tmpq]))
        mu[tmpq] <- n.trial[time] * pt
        e[tmpq] <- (y[time] - mu[tmpq])

        e.d[, tmpq] <- -n.trial[time] * pt * (1 - pt) * W.d[, tmpq]


        if (method == "NR") {

            e.dd[, , tmpq] <-  -n.trial[time] * pt * (1 - pt) * ((1 - 2 * pt) *
                              W.d[, tmpq] %o% W.d[, tmpq] + W.dd[, , tmpq])
        }
        ## update likelihood and derivatives.
        ll <- ll + y[time] * W[tmpq] - n.trial[time] * log(1 + exp(W[tmpq])) +
              log(choose(n.trial[time], y[time]))

        ll.d <- ll.d + (y[time] - mu[tmpq]) * W.d[, tmpq]
        if (method == "FS") {
            ll.dd <- ll.dd - n.trial[time] * pt * (1 - pt) * W.d[, tmpq] %o%
                     W.d[, tmpq]
        }
        if (method == "NR") {
            ll.dd <- ll.dd + (y[time] - mu[tmpq]) * W.dd[, , tmpq] -
                     n.trial[time] * pt * (1 - pt) * W.d[, tmpq] %o% W.d[, tmpq]
        }
    }
    ## end loop on time.
    list(delta = delta, ll = ll, ll.d = ll.d, ll.dd = ll.dd, eta = eta,
         W = W[mpq + 1:n], e = e[mpq + 1:n], mu = mu[mpq + 1:n],
         fitted.values = mu[mpq + 1:n])

}
