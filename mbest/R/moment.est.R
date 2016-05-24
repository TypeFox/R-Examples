# Copyright 2014 Patrick O. Perry
#
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
#
#     http://www.apache.org/licenses/LICENSE-2.0
#
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.


moment.est <- function(coefficients, nfixed, subspace, precision, dispersion,
                       start.cov = NULL, parallel = FALSE)
{
    logging::loginfo("Estimating moments", logger="mbest.mhglm.fit")
    ngroups <- nrow(coefficients)
    dim <- ncol(coefficients)
    nrandom <- dim - nfixed
    names <- colnames(coefficients)

    if (ngroups == 0L || dim == 0L){
        cov <- matrix(0, dim, dim)
        dimnames(cov) <- list(names, names)
        return(cov)
    }

    fixed <- seq_len(nfixed)
    random  <- nfixed + seq_len(nrandom)

    if (is.null(start.cov))
        start.cov <- diag(nrandom)

    logging::loginfo("Computing mean estimate and covariance bias correction",
                     logger="mbest.mhglm.fit")
    if(parallel) {
        results <- foreach(i = seq_len(ngroups)) %dopar% {
            u <- subspace[[i]]
            u1 <- u[fixed,,drop=FALSE]
            u2 <- u[random,,drop=FALSE]
            l <- precision[[i]]
            sigma2 <- dispersion[i]
            r <- length(l)
            if (r == 0L) {
                return(NULL)
            }
            s <- sqrt(l)
            us <- u %*% diag(s, r, r)
            u1s <- us[fixed,,drop=FALSE]
            u2s <- us[random,,drop=FALSE]

            cov22 <- t(u2s) %*% start.cov %*% u2s
            w.inv <- cov22 + diag(sigma2, r, r)
            R <- chol(w.inv)
            R.u1s.t <- backsolve(R, t(u1s), transpose=TRUE)
            R.u2s.t <- backsolve(R, t(u2s), transpose=TRUE)

            w11 <- t(R.u1s.t) %*% R.u1s.t
            w12 <- t(R.u1s.t) %*% R.u2s.t
            w22 <- t(R.u2s.t) %*% R.u2s.t

            R2.u2s.t <- backsolve(R, R.u2s.t)
            B <- sigma2 * t(R2.u2s.t) %*% R2.u2s.t

            w1b <- w11 %*% coefficients[i, fixed]
                    + w12 %*% coefficients[i, random]

            return(list(w11, w12, w22, w1b, B))
        }
        wtot <- (Reduce('+', lapply(results, function(x) x[[1]]))
                 / length(results))
        weight12 <- lapply(results, function(x) x[[2]])
        weight22 <- lapply(results, function(x) x[[3]])
        weight1.coef <- t(sapply(results, function(x) x[[4]]))
        wt.bias <- (Reduce('+', lapply(results, function(x) x[[5]]))
                    / length(results))
        rm(results); gc()
    } else {
        weight11 <- array(0, dim=c(ngroups, nfixed, nfixed))
        weight12 <- array(0, dim=c(ngroups, nfixed, nrandom))
        weight22 <- array(0, dim=c(ngroups, nrandom, nrandom))
        bias <- array(0, dim=c(ngroups, nrandom, nrandom))
        weight1.coef <- matrix(0, ngroups, nfixed)

        for (i in seq_len(ngroups)) {
            u <- subspace[[i]]
            u1 <- u[fixed,,drop=FALSE]
            u2 <- u[random,,drop=FALSE]
            l <- precision[[i]]
            sigma2 <- dispersion[i]
            r <- length(l)

            if (r == 0L) {
                next
            }
            # implementation trick to avoid 1/li:
            #
            # S = L^{1/2}
            #
            # W = (U2^T Sigma U2 + a L^{-1})^{-1}
            #   = [L^{-1/2} (L^{1/2} U2^T Sigma U2 L^{1/2} + a I) L^{-1/2}]^{-1}
            #   = S (S U2^T Sigma U2 S + A I)^{-1} S
            #
            # W[kl] = U[k] W U[l]^T
            #       = (U[k] S) (U2s^T Sigma U2s + a I)^{-1} (U[l] S)^T
            #
            # B = U2 W (a L^{-1}) W U2^T
            #
            #   = U2 S (S U2^T Sigma U2 S + A I)^{-1} S
            #        (a S^(-2))
            #        S (S U2^T Sigma U2 S + A I)^{-1} S U2^T
            #
            #   =  a U2 S (S U2^T Sigma U2 S + a I)^{-2} S U2^T
            #
            s <- sqrt(l)
            us <- u %*% diag(s, r, r)
            u1s <- us[fixed,,drop=FALSE]
            u2s <- us[random,,drop=FALSE]

            cov22 <- t(u2s) %*% start.cov %*% u2s
            w.inv <- cov22 + diag(sigma2, r, r)
            R <- chol(w.inv)
            R.u1s.t <- backsolve(R, t(u1s), transpose=TRUE)
            R.u2s.t <- backsolve(R, t(u2s), transpose=TRUE)

            w11 <- t(R.u1s.t) %*% R.u1s.t
            w12 <- t(R.u1s.t) %*% R.u2s.t
            w22 <- t(R.u2s.t) %*% R.u2s.t

            R2.u2s.t <- backsolve(R, R.u2s.t)
            B <- sigma2 * t(R2.u2s.t) %*% R2.u2s.t

            w1b <- (w11 %*% coefficients[i,fixed]
                    + w12 %*% coefficients[i,random])
            weight11[i,,] <- w11
            weight12[i,,] <- w12
            weight22[i,,] <- w22
            weight1.coef[i,] <- w1b
            bias[i,,] <- B
        }

        wtot <- apply(weight11, c(2,3), mean)
    }

    mean <- pseudo.solve(wtot, colMeans(weight1.coef))
    if (attr(mean, "deficient")) {
        warning(paste0("cannot solve fixed effect moment equation",
                       " due to rank deficiency"))
    }
    mean.cov <- pseudo.solve(wtot) / ngroups

    attr(mean, "deficient") <- attr(mean.cov, "deficient") <- NULL

    logging::loginfo("Computing covariance estimate", logger="mbest.mhglm.fit")
    if(parallel) {
        results <- foreach(i=seq_len(ngroups)) %dopar% {
            w12 <- weight12[[i]]
            w22 <- weight22[[i]]
            weight22.2 <- kronecker(w22, w22)
            diff <- (t(w12) %*% (coefficients[i, fixed] - mean)
                     + w22 %*% coefficients[i,random])
            weight22.coef.2 <- diff %*% t(diff)
            return(list(weight22.2, weight22.coef.2))
        }
        wtot2 <- (Reduce('+', lapply(results, function(x) x[[1]]))
                  / length(results))
        wt.cov <- (Reduce('+', lapply(results, function(x) x[[2]]))
                   / length(results))
        rm(results); gc()
    } else {
        weight22.2 <- array(NA, dim=c(ngroups, nrandom^2, nrandom^2))
        weight22.coef.2 <- array(NA, dim=c(ngroups, nrandom, nrandom))

        for (i in seq_len(ngroups)) {
            w12 <- matrix(weight12[i,,], nfixed, nrandom)
            w22 <- matrix(weight22[i,,], nrandom, nrandom)

            weight22.2[i,,] <- kronecker(w22, w22)

            diff <- (t(w12) %*% (coefficients[i,fixed] - mean)
                     + w22 %*% coefficients[i,random])
            weight22.coef.2[i,,] <- diff %*% t(diff)
        }
        wtot2 <- apply(weight22.2, c(2,3), mean)
        wt.cov <- apply(weight22.coef.2, c(2, 3), mean)
        wt.bias <- apply(bias, c(2, 3), mean)
    }

    logging::loginfo("Constructing an orthonormal basis for symmetric space",
                     logger="mbest.mhglm.fit")
    q <- nrandom
    FF <- matrix(0, q^2, q * (q + 1) / 2)
    j <- 0
    for (k in seq_len(q)) {
        for (l in seq_len(k)) {
            j <- j + 1
            f <- matrix(0, q, q)
            if (k == l) {
                f[k,l] <- 1
            } else {
                f[k,l] <- 1/sqrt(2)
                f[l,k] <- 1/sqrt(2)
            }
            FF[,j] <- as.vector(f)
        }
    }

    logging::loginfo("Solving moment equations", logger="mbest.mhglm.fit")
    tF.wtot2.F <- t(FF) %*% wtot2 %*% FF
    cov.vec <- pseudo.solve(tF.wtot2.F, t(FF) %*% as.vector(wt.cov))
    bias.vec <- pseudo.solve(tF.wtot2.F, t(FF) %*% as.vector(wt.bias))

    if (attr(cov.vec, "deficient")) {
        warning(paste0("cannot solve covariance moment equation",
                       " due to rank deficiency"))
    }

    # change back to original space
    cov <- matrix(FF %*% cov.vec, nrandom, nrandom)
    bias <- matrix(FF %*% bias.vec, nrandom, nrandom)

    # remove asymmetry arising from numerical errors
    cov <- 0.5 * (cov + t(cov))
    bias <- 0.5 * (bias + t(bias))

    eigen.cov <- eigen(cov, symmetric=TRUE)
    l <- eigen.cov$values
    u <- eigen.cov$vectors[,l > 0, drop=FALSE]
    l <- l[l > 0]
    s <- sqrt(l)
    s.u.t <- t(u) * s
    sinv.u.t <- t(u) / s
    cov.bias <- sinv.u.t %*% bias %*% t(sinv.u.t)
    eigen.cov.bias <- eigen(cov.bias, symmetric=TRUE)
    l.bias <- eigen.cov.bias$values
    u.bias.t <- t(eigen.cov.bias$vectors) %*% s.u.t
    scale <- max(1, l.bias[1])
    cov.adj <- (t(u.bias.t)
                %*% diag((scale - l.bias) / scale, length(l.bias))
                %*% u.bias.t)

    cov <- proj.psd(cov.adj)  # ensure positive definite
    if (attr(cov, "modified") || length(l) < nrow(cov) || scale != 1) {
        warning(paste("moment-based covariance matrix estimate is not positive"
                    , " semi-definite; using projection"
                    , sep=""))
    }
    attr(cov, "modified") <- NULL

    names(mean) <- names[fixed]
    dimnames(mean.cov) <- list(names[fixed], names[fixed])
    dimnames(cov) <- list(names[random], names[random])

    list(mean=mean, mean.cov=mean.cov, cov=cov)
}
