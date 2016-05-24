like.resamp <-
function(K, log.like.0, opt.cov.d, opt.mu.d, d.keep, d = 10, theta.dim = 8) 
{
    K <- K
    log.like <- NULL
    log.like.k <- log.like.0
    h.mu <- NULL
    for (i in 1:d) {
        if (!is.na(opt.cov.d[1, 1, i])) {
            h.mu <- rbind(h.mu, opt.mu.d[i, ])
        }
    }
    h.sig <- array(NA, dim = c(theta.dim, theta.dim, (K + d.keep)))
    m <- 0
    for (i in 1:d) {
        if (!is.na(opt.cov.d[1, 1, i])) {
            m <- m + 1
            h.sig[, , m] <- opt.cov.d[, , i]
        }
    }
    log.like <- log.like.k
    return(list(h.mu = h.mu, h.sig = h.sig, log.like = log.like, 
        K = K))
}

