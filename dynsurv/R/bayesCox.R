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
# Baseline prior
##############################################################################
Gamma_fun <- function(shape=0.1, rate=0.1) {
    list(shape=shape, rate=rate)
}

Const_fun <- function(value=1) {
    list(value=value, value=value)
}

GammaProcess_fun <- function(mean=0.1, ctrl=1) {
    list(mean=mean, ctrl=ctrl)
}

bp_fun <- function(type=c("Gamma", "Const", "GammaProcess"), ...) {
    type <- match.arg(type)

    if (type == "Gamma")
        hyper <- Gamma_fun(...)
    if (type == "Const")
        hyper <- Const_fun(...)
    if (type == "GammaProcess")
        hyper <- GammaProcess_fun(...)

    list(type=type, hyper=hyper)
}

##############################################################################
# Coefficient prior
##############################################################################
Normal_fun <- function(mean=0, sd=1) {
    list(mean=mean, sd=sd)
}
AR1_fun <- function(sd=1) {
    list(sd=sd, sd=sd)
}
HAR1_fun <- function(shape=2, scale=1) {
    list(shape=shape, scale=scale)
}

cp_fun <- function(type=c("Normal", "AR1", "HAR1"), ...) {
    type <- match.arg(type)

    if (type == "Normal")
        hyper <- Normal_fun(...)
    if (type == "AR1")
        hyper <- AR1_fun(...)
    if (type == "HAR1")
        hyper <- HAR1_fun(...)

    list(type=type, hyper=hyper)
}

##############################################################################
# Gibbs sampler control and general control
##############################################################################
gibbs_fun <- function(iter=3000, burn=500, thin=1, verbose=TRUE, nReport=100) {
    list(iter=iter, burn=burn, thin=thin, verbose=verbose, nReport=nReport)
}

control_bfun <- function(intercept=FALSE, a0=100, eps0=1) {
    list(intercept=intercept, a0=a0, eps0=eps0)
}

##############################################################################
# Bayesian Cox model
##############################################################################
# grid: must be sorted with last number be finite
# base.prior:
#   list(type="Gamma", shape=0.1, rate=0.1)
# coef.prior:
#   list(type="Normal", mean=0, sd=1)
#   list(type="AR1", sd=1)
#   list(type="HAR1", shape=2, scale=1)

bayesCox <- function(formula, data, grid, out,
                     model=c("TimeIndep", "TimeVarying", "Dynamic"),
                     base.prior=list(), coef.prior=list(),
                     gibbs=list(), control=list()) {

    Call <- match.call()

    # Prepare prior information
    model <- match.arg(model)
    base.prior <- do.call("bp_fun", base.prior)
    coef.prior <- do.call("cp_fun", coef.prior)

    if (model == "TimeIndep") {
        if (base.prior$type == "Gamma" && coef.prior$type == "Normal")
            id <- 11
        if (base.prior$type == "GammaProcess" && coef.prior$type == "Normal")
            id <- 12
    }
    if (model == "TimeVarying") {
        if (base.prior$type == "Gamma" && coef.prior$type == "AR1")
            id <- 21
        if (base.prior$type == "Gamma" && coef.prior$type == "HAR1")
            id <- 22
        if (base.prior$type == "GammaProcess" && coef.prior$type == "AR1")
            id <- 23
        if (base.prior$type == "GammaProcess" && coef.prior$type == "HAR1")
            id <- 24

    }
    if (model == "Dynamic") {
        if (base.prior$type == "Gamma" && coef.prior$type == "AR1")
            id <- 31
        if (base.prior$type == "Gamma" && coef.prior$type == "HAR1")
            id <- 32
        if (base.prior$type == "Const" && coef.prior$type == "AR1")
            id <- 33
        if (base.prior$type == "Const" && coef.prior$type == "HAR1")
            id <- 34
        if (base.prior$type == "GammaProcess" && coef.prior$type == "AR1")
            id <- 35
        if (base.prior$type == "GammaProcess" && coef.prior$type == "HAR1")
            id <- 36
    }

    gibbs <- do.call("gibbs_fun", gibbs)
    control <- do.call("control_bfun", control)

    # Prepare data matrix LRX
    mf <- model.frame(formula, data)
    mm <- model.matrix(formula, data)

    LRX <- cbind(mf[, 1][, 1:2], mm[, -1])
    obsInd <- which(mf[, 1][, 3] == 1)
    LRX[obsInd, 2] <- LRX[obsInd, 1]

    cov.names <- colnames(mm)[-1]

    LRX[LRX[, 2] == Inf, 2] <- max(tail(grid, 1), 999)
    LRX[is.na(LRX[, 2]), 2] <- max(tail(grid, 1), 999)

    if (control$intercept) {
        LRX <- cbind(LRX[, 1:2], 1, LRX[, -c(1:2)])
        cov.names <- c("intercept", cov.names)
    }

    colnames(LRX) <- c("L", "R", cov.names)

    # Prepare results holder
    K <- length(grid)
    nBeta <- length(cov.names)
    lambda <- rep(0, K)

    if (model == "TimeIndep")
        beta <- rep(0, nBeta)
    else
        beta <- rep(0, nBeta * K)

    nu <- rep(0, nBeta)
    jump <- rep(0, nBeta * K)

    # Call C++ function
    res <- .C("bayesCox",
              as.double(LRX), as.integer(nrow(LRX)), as.integer(ncol(LRX) - 2),
              as.double(grid), as.integer(length(grid)),
              as.character(out), as.integer(id),
              as.double(base.prior$hyper[[1]]), as.double(base.prior$hyper[[2]]),
              as.double(coef.prior$hyper[[1]]), as.double(coef.prior$hyper[[2]]),
              as.integer(gibbs$iter), as.integer(gibbs$burn), as.integer(gibbs$thin),
              as.integer(gibbs$verbose), as.integer(gibbs$nReport),
              as.double(control$a0), as.double(control$eps0),
              lambda=as.double(lambda), beta=as.double(beta), nu=as.double(nu),
              jump=as.integer(jump),
              LPML=as.double(0), DHat=as.double(0), DBar=as.double(0), pD=as.double(0),
              DIC=as.double(0))

    # Post fit processing
    if (model != "TimeIndep")
        res$beta <- matrix(res$beta, K, nBeta)

    if (coef.prior$type != "HAR1")
        res$nu <- NULL

    if (model != "Dynamic")
        res$jump <- NULL
    else
        res$jump <- matrix(res$jump, K, nBeta)

    # Return list
    rl <- list(call=Call, grid=grid, out=out, model=model, LRX=LRX,
               base.prior=base.prior, coef.prior=coef.prior,
               gibbs=gibbs, control=control,
               N=nrow(LRX), K=K, nBeta=nBeta, cov.names=cov.names,
               est=list(lambda=res$lambda, beta=res$beta, nu=res$nu, jump=res$jump),
               measure=list(LPML=res$LPML, DHat=res$DHat, DBar=res$DBar, pD=res$pD,
               DIC=res$DIC))

    class(rl) <- "bayesCox"
    rl
}
