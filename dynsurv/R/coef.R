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
# Extract the coefficient from "bayesCox" object
##############################################################################
coef.bayesCox <- function(object, ...) {

    # Monte Carlo samples
    ms <- as.matrix(read.table(file=object$out))
    dimnames(ms) <- NULL
    ms <- ms[seq(object$gibbs$burn + 1, nrow(ms), by=object$gibbs$thin), ]
    iter <- nrow(ms)

    # Dimension of baseline
    grid <- object$grid
    K <- length(grid)
    cK <- ifelse(object$model == "TimeIndep", 1, K)
    nBeta <- length(object$cov.names)

    betaMat <- as.matrix(ms[, seq(K + 1, K + nBeta * cK)])
    f <- function(x) {
        c(quantile(x, probs=0.025, names=FALSE), mean(x),
          quantile(x, probs=0.975, names=FALSE))
    }

    betaMatQT <- t(apply(betaMat, 2, f))
    if (object$model == "TimeIndep")
        betaMatQT <- betaMatQT[rep(1:nBeta, each=K), ]

    # Insert one more value at time zero
    betaMatQT <- betaMatQT[rep(seq(1, nBeta * K), rep(c(2, rep(1, K - 1)), nBeta)), ]

    res <- data.frame(betaMatQT, rep(c(0, grid), nBeta), rep(object$cov.names, each=K + 1),
                      rep(object$model, nBeta * (K + 1)))
    colnames(res) <- c("Low", "Mid", "High", "Time", "Cov", "Model")

    # Make sure the Cov retains the original order
    res$Cov <- factor(res$Cov, levels=as.character(unique(res$Cov)))

    res
}

##############################################################################
# Extract the coefficient data from "tvTran" object
##############################################################################
coef.tvTran <- function(object, ...) {
    K <- object$K
    nBeta <- object$nBeta

    rsMat <- object$rsEst[, seq(1, nBeta * K)]
    betaMat <- cbind(apply(rsMat, 2, quantile, probs=0.025, na.rm=TRUE, names=FALSE),
                     object$pEst[seq(1, nBeta * K)],
                     apply(rsMat, 2, quantile, probs=0.975, na.rm=TRUE, names=FALSE))
    # betaMat[betaMat < -bound | betaMat > bound] <- NA

    # Insert one more value at time zero
    betaMat <- betaMat[rep(seq(1, nBeta * K), rep(c(2, rep(1, K - 1)), nBeta)), ]

    res <- data.frame(betaMat, rep(c(0, object$eTime), nBeta),
                      rep(object$cov.names, each=K + 1),
                      rep("tvTran", nBeta * (K + 1)))
    colnames(res) <- c("Low", "Mid", "High", "Time", "Cov", "Model")

    # Make sure the Cov retains the original orde
    res$Cov <- factor(res$Cov, levels=as.character(unique(res$Cov)))

    res
}

##############################################################################
# Extract the coefficient from "splineCox" object
##############################################################################
coef.splineCox <- function(object, ...) {

    fit <- object$coxph.fit
    basis <- object$bsp.basis
    K <- 101

    x <- seq(basis$Boundary.knots[1], basis$Boundary.knots[2], length=K)
    bspMat <- do.call("bs", c(list(x=x), basis))

    curInd <- 1
    res <- data.frame()
    for (j in seq(1, object$nBeta)) {
        if (!object$is.tv[j]) {
            yMid <- rep(fit$coef[curInd], K)
            ySE <- sqrt(fit$var[curInd, curInd])
            curInd <- curInd + 1
        }
        else {
            sq <- seq(curInd, curInd + basis$df - 1)
            yMid <- c(bspMat %*% fit$coef[sq])
            yVar <- diag(bspMat %*% fit$var[sq, sq] %*% t(bspMat))
            yVar[which(yVar < 0)] <- 0
            ySE <- sqrt(yVar)
            curInd <- curInd + basis$df
        }

        yLow <- yMid - 1.96 * ySE
        yHigh <- yMid + 1.96 * ySE

        res <- rbind(res, data.frame(Low=yLow, Mid=yMid, High=yHigh, Time=x,
                                     Cov=object$cov.names[j], Model="Spline"))
    }

    # Make sure the Cov retains the original orde
    res$Cov <- factor(res$Cov, levels=as.character(unique(res$Cov)))

    res
}
