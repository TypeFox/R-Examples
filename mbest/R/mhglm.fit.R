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


mhglm.fit <- function(x, z, y, group, weights = rep(1, nobs),
                      start = NULL, etastart = NULL, mustart = NULL,
                      offset = rep(0, nobs), family = gaussian(),
                      control = list(), intercept = TRUE)
{
    control <- do.call("mhglm.control", control)
    x <- as.matrix(x)
    z <- as.matrix(z)
    xnames <- dimnames(x)[[2L]]
    znames <- dimnames(z)[[2L]]
    ynames <- if (is.matrix(y))
        rownames(y)
    else names(y)
    nobs <- NROW(y)
    nfixed <- ncol(x)
    nrandom <- ncol(z)
    nvars <- nfixed + nrandom
    fixed <- (seq_len(nvars) <= nfixed)
    random <- !fixed

    if (!is.null(start)) {
        if (length(start) != nfixed) {
            stop(gettextf(paste0("length of 'start' should equal %d",
                                 " and correspond to initial coefs for %s"),
                          nfixed, paste(deparse(xnames), collapse=", ")),
                 domain=NA)
        }
        start <- c(start, numeric(nrandom))
    }

    # group-specific estimates
    m <- rdglm.group.fit(x = cbind(x, z), y = y, group = group,
                         weights = weights, start = start,
                         etastart = etastart, mustart = mustart,
                         offset = offset, family = family, 
                         parallel = control$parallel,
                         control = control$fit.control,
                         method = control$fit.method,
                         intercept = intercept)
    ngroups <- m$ngroups

    # compute pooled dispersion estimates
    dispersion.tot <- dispersion.pooled(m$dispersion, m$df.residual)
    df.residual.tot <- sum(m$df.residual)
    dispersion <- rep(dispersion.tot, ngroups)

    if (dispersion.tot == 0)
        stop("cannot estimate dispersion (no residual degrees of freedom)")

    # compute group-sqecific precision square root (unpivoted)
    Rp <- as.list(rep(NULL, ngroups))

    if(control$parallel) {
        Rp <- foreach(i = seq_len(ngroups)) %dopar% {
            qr.i <- m$qr[[i]]
            rank.i <- qr.i$rank
            Rp <- matrix(0, rank.i, nvars)
            if (rank.i > 0L) {
                R.i <- qr.R(qr.i)
                pivot.i <- qr.i$pivot
                Rp[, pivot.i] <- R.i[seq_len(rank.i), ]
            }
            return(Rp)
        }
    } else {
        for(i in seq_len(ngroups)) {
            qr.i <- m$qr[[i]]
            rank.i <- qr.i$rank
            Rp[[i]] <- matrix(0, rank.i, nvars)
            if (rank.i > 0L) {
                R.i <- qr.R(qr.i)
                pivot.i <- qr.i$pivot
                Rp[[i]][,pivot.i] <- R.i[seq_len(rank.i),]
            }
        }
    }

    # change coordinates so that average precision is identity

    if (control$standardize) {
        # compute averge precision of estimates
        prec.avg <- matrix(0, nvars, nvars)
        if(control$parallel) {
            prec.avg <- foreach(i = seq_len(ngroups)) %dopar% {
                scale.i <- sqrt(1/dispersion[i])
                return(crossprod(scale.i * Rp[[i]]))
            }
            prec.avg <- Reduce('+', prec.avg)
        } else {
            for (i in seq_len(ngroups)) {
                scale.i <- sqrt(1/dispersion[i])
                prec.avg <- prec.avg + crossprod(scale.i * Rp[[i]])
            }
        }
        prec.avg <- prec.avg / ngroups

        # cov(coef) = (t(R) R)^{-1} = R^{-1} R^{-T}
        # cov(R x) = R cov(x) R^T
        # [cov(R x)]^{-1} = R^{-T} [cov(x)]^{-1} R^{-1}
        suppressWarnings({
            R.fixed <- chol(prec.avg[fixed,fixed], pivot = TRUE)
            R.random <- chol(prec.avg[random,random], pivot = TRUE)
        })

        #pivot <- attr(R, "pivot")
        #rank <- attr(R, "rank")
        #r1.fixed <- seq_len(attr(R.fixed, "rank"))
        #r1.random <- seq_len(attr(R.random, "rank"))
        #R.fixed <- R.fixed[r1.fixed,r1.fixed,drop=FALSE]
        #R.random <- R.random[r1.random,r1.random,drop=FALSE]
    }
    else { # control.standardize==FALSE
        R.fixed <- diag(nfixed)
        attr(R.fixed, "pivot") <- seq_len(nfixed)
        attr(R.fixed, "rank") <- nfixed
        R.random <- diag(nrandom)
        attr(R.random, "pivot") <- seq_len(nrandom)
        attr(R.random, "rank") <- nrandom

        #pivot <- seq_len(nvars)
        #rank <- nvars
        #r1 <- seq_len(rank)
    }

    rank.fixed <- attr(R.fixed, "rank")
    rank.random <- attr(R.random, "rank")
    rank <- rank.fixed + rank.random

    r1.fixed <- seq_len(rank.fixed)
    r1.random <- seq_len(rank.random)
    r1 <- seq_len(rank)

    pivot.fixed <- which(fixed)[attr(R.fixed, "pivot")]
    pivot.random <- which(random)[attr(R.random, "pivot")]
    pivot <- c(pivot.fixed[r1.fixed],  pivot.random[r1.random],
               pivot.fixed[-r1.fixed], pivot.random[-r1.random])

    R <- matrix(0, rank, rank)
    R[r1.fixed, r1.fixed] <- R.fixed[r1.fixed, r1.fixed]
    (R[rank.fixed + r1.random, rank.fixed + r1.random]
        <- R.random[r1.random, r1.random])

    # compute standardized coeficients:
    #   coef[i,] <- R %*% m$coefficients[i,pivot[r1]]
    coef <- m$coefficients[,pivot[r1],drop=FALSE] %*% t(R)

    # compute group-specific subspace and precisions
    if(control$parallel) {
        results <- foreach(i = seq_len(ngroups)) %dopar% {
            if (nrow(Rp[[i]]) > 0L) {
                prec.sqrt <- backsolve(R, t(Rp[[i]][,pivot[r1],drop=FALSE]),
                                       transpose=TRUE)
                prec.sqrt.svd <- svd(prec.sqrt)
                return(list(subspace=prec.sqrt.svd$u,
                            precision=(prec.sqrt.svd$d)^2))
            } else {
                return(list(subspace=matrix(0, rank, 0), precision=numeric()))
            }
        }
        subspace <- lapply(results, function(x) x[['subspace']])
        precision <- lapply(results, function(x) x[['precision']])
        rm(results)
    } else {
        subspace <- as.list(rep(NULL, ngroups))
        precision <- as.list(rep(NULL, ngroups))
        for (i in seq_len(ngroups)) {
            if (nrow(Rp[[i]]) > 0L) {
                prec.sqrt <- backsolve(R, t(Rp[[i]][,pivot[r1],drop=FALSE]),
                                       transpose=TRUE)
                prec.sqrt.svd <- svd(prec.sqrt)
                subspace[[i]] <- prec.sqrt.svd$u
                precision[[i]] <- (prec.sqrt.svd$d)^2
            } else {
                subspace[[i]] <- matrix(0, rank, 0)
                precision[[i]] <- numeric()
            }
        }
    }

    # compute coefficient mean and covariance estimates
    logging::loginfo("Computing mean and covariance estimates",
                     logger="mbest.mhglm.fit")
    est0 <- NULL
    suppressWarnings({
        for (s in seq_len(control$steps)) {
            est0 <- moment.est(coef, nfixed=rank.fixed, subspace, precision,
                               dispersion, start.cov=est0$cov, 
                               parallel=control$parallel)
            logging::loginfo("Refining mean and covariance estimates",
                             logger="mbest.mhglm.fit")
        }
    })
    est <- moment.est(coef, nfixed=rank.fixed, subspace, precision, dispersion,
                      start.cov=est0$cov, parallel=control$parallel)
    mean <- est$mean
    mean.cov <- est$mean.cov
    cov <- est$cov

    # change back to original coordinates
    coef.mean <- rep(NA, nfixed)
    coef.mean.cov <- matrix(NA, nfixed, nfixed)
    coef.cov <- matrix(NA, nrandom, nrandom)
    (coef.mean[attr(R.fixed, "pivot")[r1.fixed]]
        <- backsolve(R.fixed, mean))
    (coef.mean.cov[attr(R.fixed, "pivot")[r1.fixed],
                   attr(R.fixed, "pivot")[r1.fixed]]
        <- backsolve(R.fixed, t(backsolve(R.fixed, mean.cov))))
    (coef.cov[attr(R.random, "pivot")[r1.random],
              attr(R.random, "pivot")[r1.random]]
        <- backsolve(R.random, t(backsolve(R.random, cov))))

    # set coordinate names
    names(coef.mean) <- xnames
    dimnames(coef.mean.cov) <- list(xnames, xnames)
    dimnames(coef.cov) <- list(znames, znames)

    fit <- list(family = family, coefficient.mean = coef.mean,
                coefficient.mean.cov = coef.mean.cov,
                coefficient.cov = coef.cov, coefficients = coef,
                subspace = subspace, precision = precision,
                dispersion = dispersion.tot, df.residual = df.residual.tot,
                R = R, rank = rank, rank.fixed = rank.fixed,
                rank.random = rank.random, pivot = pivot, y = y, group = group,
                prior.weights = weights, offset = offset, nobs = nobs)
    fit
}

