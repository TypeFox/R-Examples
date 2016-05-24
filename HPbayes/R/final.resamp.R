final.resamp <- function (K, B1, H.new, H.k, log.like, d.keep, prior, h.mu, h.sig, 
    nrisk, ndeath, B = 400, theta.dim = 8, age = c(1e-05, 1, 
        seq(5, 100, 5))) 
{
    lx <- nrisk
    dx <- ndeath
    mp8.mle <- function(theta, x.fit = age) {
        p.hat <- mod8p(theta = theta, x = x.fit)
        ll <- ll.binom(x = dx, n = lx, p = p.hat)
        return(ll)
    }
    B0 <- theta.dim * 1000
    prior.cov <- cov(prior)
    imp.keep <- theta.dim * 100
    vwts <- NULL
    ewts <- NULL
    mwts <- NULL
    nup <- NULL
    nck <- NULL
    mwt.case <- NULL
    frac.nup <- NULL
    prilo <- apply(prior, 2, min)
	prihi <- apply(prior, 2, max)
	lo <- prilo
	hi <- prihi

    for (k in 1:K) {
        if (k == 1) {
            log.like.k <- rep(NA, B1)
            for (i in 1:B1) {
                log.like.k[i] <- mp8.mle(H.new[i, ])
                if (i%%100 == 0) {
                }
            }
        }
        if (k > 1) {
            log.like.k <- rep(NA, B)
            for (i in 1:B) {
                log.like.k[i] <- mp8.mle(H.new[i, ])
                if (i%%100 == 0) {
                }
            }
            H.k <- rbind(H.k, H.new)
        }
        log.like <- c(log.like, log.like.k)
        n.k <- B * (k - 1) + B0 + B1
        ncol.qk <- 1 + d.keep + (k - 1)
        DENS.MAT <- matrix(NA, nrow = n.k, ncol = ncol.qk)
        DENS.MAT[, 1] <- (B0/n.k) * dens.prior(H.k, pri.lo = prilo, 
            pri.hi = prihi)
        for (i in 1:(ncol.qk - 1)) {
            DENS.MAT[, (i + 1)] <- B/n.k * dmvnorm(x = H.k, mean = h.mu[i, 
                ], sigma = h.sig[, , i])
        }
        qk <- rowSums(DENS.MAT)
        a <- -max(log.like, na.rm = TRUE)
        like.k <- exp(log.like)
        wts.k.0 <- like.k * dens.prior(H.k, pri.lo = prilo, pri.hi = prihi)
        wts.k <- wts.k.0/qk
        wts.k <- wts.k/sum(wts.k, na.rm = TRUE)
        wts.k[is.na(wts.k)] <- 0
        maxw.k <- max(wts.k)
        thresh.k <- min(maxw.k, 0.05)
        imp.theta <- H.k[wts.k == maxw.k, ]
        dist.to.mean <- mahalanobis(x = H.k, center = imp.theta, 
            cov = prior.cov)
        keep.4.cov <- sort(dist.to.mean, decreasing = FALSE)[imp.keep]
        inputs.4.cov <- H.k[dist.to.mean <= keep.4.cov, ]
        wts.4.cov <- wts.k[dist.to.mean <= keep.4.cov]
        wts.4.cov <- ((wts.4.cov + 1/length(wts.k)) * 0.5)/sum(((wts.4.cov + 
            1/length(wts.k)) * 0.5))
        Sigma.k <- cov.wt(inputs.4.cov, wt = wts.4.cov)$cov
        h.mu <- rbind(h.mu, imp.theta)
        h.sig[, , (d.keep + k)] <- Sigma.k
        H.new <- mvrnorm(n = B, mu = imp.theta, Sigma = Sigma.k)
        for (i in 1:B) {
            lt0 <- table(H.new[i, ] < 0)
            while (lt0[1] != ncol(prior) 
            | H.new[i, 1] < lo[1] | H.new[i, 1] > hi[1] 
            | H.new[i, 2] < lo[2] | H.new[i, 2] > hi[2] 
            | H.new[i, 3] < lo[3] | H.new[i, 3] > hi[3] 
            | H.new[i, 4] < lo[4] | H.new[i, 4] > hi[4] 
            | H.new[i, 5] < lo[5] | H.new[i, 5] > hi[5] 
            | H.new[i, 6] < lo[6] | H.new[i, 6] > hi[6] 
            | H.new[i, 7] < lo[7] | H.new[i, 7] > hi[7] 
            | H.new[i, 8] < lo[8] | H.new[i, 8] > hi[8]) {
                H.new[i, ] <- mvrnorm(n = 1, mu = imp.theta, 
                  Sigma = Sigma.k)
                lt0 <- table(H.new[i, ] < 0)
            }
            if (i%%100 == 0) {
            }
        }
        H.new <- cbind(H.new[, 1:8])
        vwts <- c(vwts, var.rwts(wts.k))
        ewts <- c(ewts, entropy.wts(wts.k))
        mwts <- c(mwts, maxw.k)
        nup <- c(nup, expt.upts(wts.k, m = B))
        frac.nup <- c(frac.nup, expt.upts(wts.k, m = B)/B)
        N.k <- K * B
        case.n <- 1:N.k
        mwt.case <- c(mwt.case, case.n[wts.k == maxw.k])
        assign(paste("imp.theta.k", k, sep = ""), imp.theta)
        assign(paste("Sigma.k", k, sep = ""), Sigma.k)
        if (k%%5 == 0) {
            assign(paste("log.like.k", k, sep = ""), log.like.k)
            assign(paste("wts.k", k, sep = ""), wts.k)
            assign(paste("maxw.k", k, sep = ""), maxw.k)
            assign(paste("imp.theta", k, sep = ""), imp.theta)
            assign(paste("H.k", k, sep = ""), H.k)
            assign(paste("qk", k, sep = ""), qk)
            save.vectors <- c(paste("log.like.k", k, sep = ""), 
                paste("wts.k", k, sep = ""), paste("maxw.k", 
                  k, sep = ""), paste("imp.theta", k, sep = ""), 
                paste("Sigma.k", k, sep = ""), paste("H.k", k, 
                  sep = ""), "qk", "log.like")
            rm.vectors <- c(paste("log.like.k", k, sep = ""), 
                paste("wts.k", k, sep = ""), paste("maxw.k", 
                  k, sep = ""))
            rm(list = rm.vectors)
        }
        if (frac.nup[k] >= 0.632) {
            break
        }
    }
    return(list(H.new = H.new, vwts = vwts, ewts = ewts, mwts = mwts, 
        nup = nup, frac.nup = frac.nup, wts.k = wts.k, mwt.case=mwt.case))
}