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


ebayes.est <- function(coefficients, nfixed, subspace, precision, dispersion,
                       coefficient.mean, coefficient.cov,
                       coefficient.cov.sqrt, postVar = FALSE)
{
    coef <- coefficients
    coef.mu <- coefficient.mean
    coef.cov <- coefficient.cov
    r <- length(precision)

    if (postVar) {
        R <- coefficient.cov.sqrt
        pivot <- attr(R, "pivot")
        rank <- attr(R, "rank")
        R <- R[seq_len(rank),,drop=FALSE]
    }

    nrandom <- nrow(coef.cov)
    fixed <- seq_len(nfixed)
    random <- nfixed + seq_len(nrandom)

    if (r == 0L) {
        coef.eb <- numeric(nrandom)
        cov.eb <- coef.cov
    } else {
        # implementation trick to avoid 1/li:
        # U (U^T Sigma U + a L^{-1})^{-1} U^T
        #   = U L^{1/2} (L^{1/2} U^T Sigma U L^{1/2} + a I)^{-1} L^{1/2} U^T
        #   = Us (Us^T Sigma Us + a I)^{-1} Us^T
        u <- subspace
        s <- sqrt(precision)
        us <- u %*% diag(s, r, r)
        u1s <- us[fixed,,drop=FALSE]
        u2s <- us[random,,drop=FALSE]

        cov.ii <- t(u2s) %*% coef.cov %*% u2s
        h <- cov.ii + diag(dispersion, r, r)
        w.diff <- u2s %*% solve(h, (t(u1s) %*% (coef[fixed] - coef.mu)
                                    + t(u2s) %*% coef[random]))
        coef.eb <- coef.cov %*% w.diff

        if (postVar) {
            Ru2s <- R %*% u2s[pivot,,drop=FALSE]
            C <- chol(diag(dispersion, rank, rank) + tcrossprod(Ru2s))
            R.post <- matrix(0, rank, nrandom)
            R.post[,pivot] <- backsolve(C, R, transpose=TRUE)
            cov.eb <- dispersion * crossprod(R.post)
        }
    }

    if (postVar)
        attr(coef.eb, "postVar") <- cov.eb

    coef.eb
}


ebayes.group.est <- function(coefficients, nfixed, subspace, precision, dispersion,
                             coefficient.mean, coefficient.cov, postVar = FALSE)
{
    ngroups <- nrow(coefficients)
    nvars <- ncol(coefficients) - nfixed

    coefficients.eb <- matrix(NA, ngroups, nvars)
    dimnames(coefficients.eb) <- dimnames(coefficients)

    if (postVar) {
        suppressWarnings({
            coefficient.cov.sqrt <- chol(coefficient.cov, pivot=TRUE)
        })
    } else {
        coefficient.cov.sqrt <- NULL
    }

    if (postVar) {
        cov.eb <- array(NA, c(nvars, nvars, ngroups))
        dimnames(cov.eb) <- list(colnames(coefficients),
                                 colnames(coefficients),
                                 rownames(coefficients))
    }

    for (i in seq_len(ngroups)) {
        eb <- ebayes.est(coefficients[i,],
                         nfixed = nfixed,
                         subspace = subspace[[i]],
                         precision = precision[[i]], dispersion = dispersion[i],
                         coefficient.mean = coefficient.mean,
                         coefficient.cov = coefficient.cov,
                         coefficient.cov.sqrt = coefficient.cov.sqrt,
                         postVar = postVar)

        coefficients.eb[i,] <- eb
        if (postVar) {
            cov.eb[,,i] <- attr(eb, "postVar")
        }
    }

    if (postVar)
        attr(coefficients.eb, "postVar") <- cov.eb

    coefficients.eb
}

