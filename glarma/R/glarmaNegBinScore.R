glarmaNegBinScore <- function(y, X, offset = NULL, delta, phiLags, thetaLags,
                              method  = "FS") {
    r <- ncol(X)
    n <- length(y)
    p <- length(phiLags)
    q <- length(thetaLags)
    s <- r + p + q + 1
    beta <- delta[1:r]
    phi <- delta[(r + 1):(r + p)]
    theta <- delta[(r + p + 1):(r + p + q)]
    alpha <- delta[s]
    u <- c(rep(0, s - 1), 1)
    mpq <- 0

    if ((p + q) > 0) {
        mpq <- max(phiLags[p], thetaLags[q])
    }

    nmpq <- n + mpq
    e <- array(0, nmpq)
    Z <- array(0, nmpq)
    W <- array(0, nmpq)
    mu <- array(0, nmpq)
    e.d <- array(0, c(s, nmpq))
    Z.d <- array(0, c(s, nmpq))
    W.d <- array(0, c(s, nmpq))
    if (method  == "NR") {
        e.dd <- array(0, c(s, s, nmpq))
        Z.dd <- array(0, c(s, s, nmpq))
        W.dd <- array(0, c(s, s, nmpq))
    }
    if(is.null(offset)) eta<-X %*% beta else eta<- X %*% beta + offset
    ll <- 0
    ll.d <- matrix(0, ncol = 1, nrow = s)
    ll.dd <- matrix(0, ncol = s, nrow = s)

    for (time in 1:n) {
        tmpq <- time + mpq

        if (p > 0) {
            Z.d[(r + 1):(r + p), tmpq] <- Z[tmpq - phiLags] + e[tmpq - phiLags]
            if (method  == "NR") {
                Z.dd[(r + 1):(r + p), , tmpq] <-
                    t((Z.d + e.d)[, (tmpq - phiLags)])
                Z.dd[, (r + 1):(r + p), tmpq] <-
                    Z.dd[, (r + 1):(r + p), tmpq] +
                        (Z.d + e.d)[, (tmpq - phiLags)]
            }
            for (i in 1:p) {
                Z[tmpq] <- Z[tmpq] + phi[i] * (Z + e)[tmpq - phiLags[i]]
                Z.d[, tmpq] <- Z.d[, tmpq] + phi[i] *
                               (Z.d[, tmpq - phiLags[i]] +
                               e.d[, tmpq - phiLags[i]])
                if (method  == "NR") {
                  Z.dd[, , tmpq] <- Z.dd[, , tmpq] + phi[i] *
                                    (Z.dd[, , tmpq - phiLags[i]] +
                                    e.dd[, , tmpq - phiLags[i]])
                }
            }
        }

        if (q > 0) {

            Z.d[(r + p + 1):(r + p + q), tmpq] <- e[tmpq - thetaLags]
            if (method  == "NR") {
                Z.dd[(r + p + 1):(r + p + q), , tmpq] <-
                    Z.dd[(r + p + 1):(r + p + q), , tmpq] +
                        t(e.d[, tmpq - thetaLags])
                Z.dd[, (r + p + 1):(r + p + q), tmpq] <-
                    Z.dd[, (r + p + 1):(r + p + q), tmpq] +
                        e.d[, tmpq - thetaLags]
            }
            for (i in 1:q) {

                Z[tmpq] <- Z[tmpq] + theta[i] * e[tmpq - thetaLags[i]]
                Z.d[, tmpq] <- Z.d[, tmpq] + theta[i] *
                               e.d[, tmpq - thetaLags[i]]
                if (method  == "NR") {
                  Z.dd[, , tmpq] <- Z.dd[, , tmpq] + theta[i] *
                                    e.dd[, , tmpq - thetaLags[i]]
                }
            }
        }

        W[tmpq] <- eta[time] + Z[tmpq]
        W.d[, tmpq] <- matrix(c(X[time, ], rep(0, p + q + 1)), ncol = 1) +
                       Z.d[, tmpq]
        if (method  == "NR") {
            W.dd[, , tmpq] <- Z.dd[, , tmpq]
        }
        mu[tmpq] <- exp(W[tmpq])

        mut <- mu[tmpq]
        yt <- y[time]
        var <- mut + mut^2/alpha


        e.W <- (-mut * (alpha * yt + mut * (alpha + 2 * yt)))/
               (2 * alpha * var^(3/2))
        e.a <- mut^2 * (yt - mut)/(2 * alpha^2 * var^(3/2))
        if (method  == "NR") {
            e.WW <- (yt * alpha^2 + mut * alpha * (2 * yt - alpha) +
                     2 * mut^2 *
                    (2 * yt + alpha))/(4 * (mut + alpha)^2 * var^(1/2))
            e.aW <- mut^3 * (alpha * yt - mut * (2 * yt + 3 * alpha))/
                    (4 * alpha^3 * var^(5/2))
            e.aa <- mut^3 * (mut - yt) * (mut + 4 * alpha)/
                    (4 * alpha^4 * var^(5/2))
        }
        e[tmpq] <- (yt - mut)/mut
        e.d[, tmpq] <- -(1 + e[tmpq]) * W.d[, tmpq]
        if (method  == "NR") {
            e.dd[, , tmpq] <- -(1 + e[tmpq]) * W.dd[, , tmpq] -
                               (e.d[, tmpq]) %o% W.d[, tmpq]
        }

        ## update likelihood and derivatives.

        ll <- ll + log(gamma(alpha + yt)/(gamma(alpha) * gamma(yt + 1))) +
              alpha * log(alpha/(alpha + mut)) + yt * log(mut/(alpha + mut))
        ll.W <- alpha * (yt - mut)/(mut + alpha)
        ll.a <- (digamma(alpha + yt) - digamma(alpha) + (mut - yt)/
                (mut + alpha) + log(alpha/(mut + alpha)))
        ll.d <- ll.d + ll.W * W.d[, tmpq] + ll.a * u
        if (method  == "FS") {
            ll.dd <- ll.dd + (-(alpha * mut * (yt + alpha)/(mut + alpha)^2) *
                     W.d[, tmpq] %o% W.d[, tmpq] + (trigamma(alpha + yt) -
                     trigamma(alpha) + (yt * alpha + mut^2)/
                     (alpha * (alpha + mut)^2)) * u %o% u)
        }

        if (method  == "NR") {
            ll.dd <- ll.dd + (alpha * (yt - mut)/(mut + alpha) *
                     W.dd[, , tmpq] - (alpha * mut * (yt + alpha)/
                     (mut + alpha)^2) * W.d[, tmpq] %o% W.d[, tmpq] +
                     ((yt - mut) * mut/(mut + alpha)^2) * (W.d[, tmpq] %o% u +
                     u %o% W.d[, tmpq]) + (trigamma(alpha + yt) -
                     trigamma(alpha) + (yt * alpha + mut^2)/
                     (alpha * (alpha + mut)^2)) * u %o% u)
        }
    }
    list(delta = delta, ll = ll, ll.d = ll.d, ll.dd = ll.dd, eta = eta,
         W = W[mpq + 1:n], e = e[mpq + 1:n], mu = mu[mpq + 1:n],
         fitted.values = mu[mpq + 1:n])

}
