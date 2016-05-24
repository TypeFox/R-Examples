em.hmm <-
function (x, L = 2, maxiter = 1000, est.null = FALSE) 
{
    NUM <- length(x)
    ptol <- 1e-04
    niter <- 0
    converged = TRUE
    if (L == 1) {
        pii.new <- c(0.5, 0.5)
        A.new <- matrix(c(0.8, 0.2, 0.4, 0.6), 2, 2, byrow = TRUE)
        f0.new <- c(2, 1)
        f1.new <- c(4, 1)
        diff <- 1
        while (diff > ptol && niter < maxiter) {
            niter <- niter + 1
            pii.old <- pii.new
            A.old <- A.new
            f0.old <- f0.new
            f1.old <- f1.new
            bwfw.res <- bwfw1.hmm(x, pii.old, A.old, f0.old, 
                f1.old)
            gamma <- bwfw.res$pr
            dgamma <- bwfw.res$ts
            c0 <- bwfw.res$rescale
            for (i in 0:1) {
                pii.new[i + 1] <- gamma[1, i + 1]
            }
            for (i in 0:1) {
                for (j in 0:1) {
                  q1 <- sum(dgamma[i + 1, j + 1, ])
                  q2 <- sum(gamma[1:(NUM - 1), i + 1])
                  A.new[i + 1, j + 1] <- q1/q2
                }
            }
            q5 <- sum(gamma[, 1] * x)
            mu0 <- q5/sum(gamma[, 1])
            q6 <- sum(gamma[, 1] * (x - mu0) * (x - mu0))
            sd0 <- sqrt(q6/sum(gamma[, 1]))
            f0.new <- c(mu0, sd0)
            if (!est.null) {
                f0.new <- c(0, 1)
            }
            q1 <- sum(gamma[, 2])
            q2 <- sum(gamma[, 2] * x)
            mu1 <- q2/q1
            q3 <- sum(gamma[, 2] * (x - mu1) * (x - mu1))
            sd1 <- sqrt(q3/q1)
            f1.new <- c(mu1, sd1)
            df1 <- abs(A.old - A.new)
            df2 <- abs(f1.old - f1.new)
            diff <- max(df1, df2)
            if (is.na(diff)) {
                converged = FALSE
                break
            }
        }
        lfdr <- gamma[, 1]
        if (converged) {
            logL <- -sum(log(c0))
            if (est.null) {
                BIC <- logL - (3 * L + 4) * log(NUM)/2
            }
            else {
                BIC <- logL - (3 * L + 2) * log(NUM)/2
            }
            em.var <- list(pii = pii.new, A = A.new, f0 = f0.new, 
                f1 = f1.new, LIS = lfdr, logL = logL, BIC = BIC, 
                ni = niter, converged = converged)
        }
        else {
            BIC <- logL <- (-Inf)
            em.var <- list(pii = pii.old, A = A.old, f0 = f0.old, 
                f1 = f1.old, LIS = lfdr, logL = logL, BIC = BIC, 
                ni = niter, converged = converged)
        }
    }
    else if (L > 1) {
        pii.new <- c(0.5, 0.5)
        A.new <- matrix(c(0.95, 0.05, 0.5, 0.5), 2, 2, byrow = TRUE)
        pc.new <- rep(1, L)/L
        mus <- seq(from = -1, by = 1.5, length = L)
        sds <- rep(1, L)
        f0.new <- c(2, 1)
        f1.new <- cbind(mus, sds)
        diff <- 1
        while (diff > ptol && niter < maxiter) {
            niter <- niter + 1
            pii.old <- pii.new
            A.old <- A.new
            pc.old <- pc.new
            f0.old <- f0.new
            f1.old <- f1.new
            bwfw.res <- bwfw.hmm(x, pii.old, A.old, pc.old, f0.old, 
                f1.old)
            gamma <- bwfw.res$pr
            dgamma <- bwfw.res$ts
            omega <- bwfw.res$wt
            c0 <- bwfw.res$rescale
            for (i in 0:1) {
                pii.new[i + 1] <- gamma[1, i + 1]
            }
            for (i in 0:1) {
                for (j in 0:1) {
                  q1 <- sum(dgamma[i + 1, j + 1, ])
                  q2 <- sum(gamma[1:(NUM - 1), i + 1])
                  A.new[i + 1, j + 1] <- q1/q2
                }
            }
            q5 <- sum(gamma[, 1] * x)
            mu0 <- q5/sum(gamma[, 1])
            q6 <- sum(gamma[, 1] * (x - mu0) * (x - mu0))
            sd0 <- sqrt(q6/sum(gamma[, 1]))
            f0.new <- c(mu0, sd0)
            if (!est.null) {
                f0.new <- c(0, 1)
            }
            mus <- 1:L
            sds <- 1:L
            for (c in 1:L) {
                q1 <- sum(omega[, c])
                q2 <- sum(gamma[, 2])
                pc.new[c] <- q1/q2
                q3 <- sum(omega[, c] * x)
                mus[c] <- q3/q1
                q4 <- sum(omega[, c] * (x - mus[c]) * (x - mus[c]))
                sds[c] <- sqrt(q4/q1)
            }
            f1.new <- cbind(mus, sds)
            df1 <- abs(A.old - A.new)
            df2 <- abs(f1.old - f1.new)
            diff <- max(df1, df2)
            if (is.na(diff)) {
                converged = FALSE
                break
            }
        }
        lfdr <- gamma[, 1]
        if (converged) {
            logL <- -sum(log(c0))
            if (est.null) {
                BIC <- logL - (3 * L + 4) * log(NUM)/2
            }
            else {
                BIC <- logL - (3 * L + 2) * log(NUM)/2
            }
            em.var <- list(pii = pii.new, A = A.new, pc = pc.new, 
                f0 = f0.new, f1 = f1.new, LIS = lfdr, logL = logL, 
                BIC = BIC, ni = niter, converged = converged)
        }
        else {
            logL <- (-Inf)
            BIC <- logL <- (-Inf)
            em.var <- list(pii = pii.old, A = A.old, pc = pc.old, 
                f0 = f0.old, f1 = f1.old, LIS = lfdr, logL = logL, 
                BIC = BIC, ni = niter, converged = converged)
        }
    }
    return(em.var)
}
