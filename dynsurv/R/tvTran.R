##############################################################################
##
##   R package dynsurv by Xiaojing Wang, Jun Yan, and Ming-Hui Chen
##   Copyright (C) 2011
##
##   This file is part of the R package dynsurv.
##
##   The R package dynsurv is free software: you can redistribute it and/or
##   modify it under the terms of the GNU General Public License as published
##   by the Free Software Foundation, either version 3 of the License, or
##   (at your option) any later version.
##
##   The R package dynsurv is distributed in the hope that it will be useful,
##   but WITHOUT ANY WARRANTY; without even the implied warranty of
##   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
##   GNU General Public License for more details.
##
##   You should have received a copy of the GNU General Public License
##   along with the R package dynsurv. If not, see <http://www.gnu.org/licenses/>.
##
##############################################################################

##############################################################################
# Utility functions
##############################################################################
control_tfun <- function(resample=TRUE, R=30) {
    list(resample=resample, R=R)
}

tvTran_lite <- function(X, dNMat, YMat, offSetMat) {

    # Estimating function
    f <- function(beta, preExpXb, dN, Y, offSet) {
        c(t(X) %*% (dN - Y * (exp(c(X %*% beta)) - preExpXb))) + offSet
    }

    N <- nrow(X)
    nBeta <- ncol(X)
    K <- ncol(dNMat)

    # Initialization
    b0 <- rep(0, nBeta)
    betaMat <- matrix(0, K, nBeta)
    preExpXb <- rep(0, N)
    termCode <- rep(0, K)

    # Solve estimating equations sequentially
    for (k in 1:K) {
        dNVec <- dNMat[, k]
        YVec <- YMat[, k]

        # Solve time-varying coefficients at time t_i
        res <- nleqslv(b0, f, preExpXb=preExpXb, dN=dNMat[, k], Y=YMat[, k],
                       offSet=offSetMat[k, ], xscalm="auto")
        beta <- res$x
        termCode[k] <- res$termcd
        betaMat[k, ] <- beta
        preExpXb <- exp(c(X %*% beta))
    }

    betaMat[termCode != 1, ] <- NA

    # Append termination codes to the end of coefficient estimates
    c(c(betaMat), termCode)
}

##############################################################################
# Time-varying coefficient transformation model, Peng and Huang (2007)
##############################################################################
tvTran <- function(formula, data, control=list()) {

    Call <- match.call()
    control <- do.call("control_tfun", control)

    # Right censored data
    mf <- model.frame(formula, data)
    rsp <- mf[, 1]

    # Event subject index
    eIndex <- which(rsp[, "status"] == 1)

    # Unique event time
    eTime <- sort(unique(rsp[eIndex, "time"]))
    K <- length(eTime)

    # Design matrix with intercept term
    X <- model.matrix(formula, data)
    N <- nrow(X)
    nBeta <- ncol(X)
    cov.names <- c("intercept", colnames(X)[-1])
    X <- matrix(X, N, nBeta)

    # Prepare event, at-risk and offset matrix
    dNMat <- matrix(0, N, K)
    dNMat[eIndex, ] <- outer(rsp[eIndex, "time"], eTime, "==") + 0

    YMat <- outer(rsp[, "time"], eTime, ">=") + 0

    offSetMat <- matrix(0, K, nBeta)

    # Point estimate
    pEst <- tvTran_lite(X, dNMat, YMat, offSetMat)

    # Resampling
    if (control$resample) {
        R <- control$R
        rsEst <- matrix(0, R, (nBeta + 1) * K)

        # Indicator matrix: I(time_i <= eTime_j) * I(status_i == 1)
        indMat <- outer(rsp[, "time"], eTime, "<=") * rsp[, "status"]

        for (r in 1:R) {
            zeta <- rnorm(N)
            offSetMat <- apply(cbind(0, t(X) %*% (indMat * zeta)), 1, diff)
            rsEst[r, ] <- tvTran_lite(X, dNMat, YMat, offSetMat)
        }
    }
    else
        reEst <- NULL

    # Return list
    rl <- list(call=Call, eTime=eTime, control=control,
               N=N, K=K, nBeta=nBeta, cov.names=cov.names,
               pEst=pEst, rsEst=rsEst)
    class(rl) <- "tvTran"

    rl
}
