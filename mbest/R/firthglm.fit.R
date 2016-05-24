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


# improve.tol = 1e-4 suggested by Nocedal & Wright (c1, p. 38)
# curvature.tol = 0.9 suggested by Nocedal & Wright (c2, p. 39)

firthglm.control <- function(epsilon = 1e-8, maxit = 25, qr.tol = 1e-7,
                             improve.tol = 1e-4, curvature.tol = 0.9,
                             linesearch.method = "linesearch",
                             linesearch.maxit = 20, trace = FALSE)
{
    if (!is.numeric(epsilon) || epsilon <= 0)
        stop("value of 'epsilon' must be > 0")
    if (!is.numeric(maxit) || maxit <= 0)
        stop("maximum number of iterations must be > 0")
    if (!is.numeric(qr.tol) || qr.tol <= 0)
        stop("value of 'qr.tol' must be > 0")
    if (!is.numeric(improve.tol) || improve.tol <= 0)
        stop("value of 'improve.tol' must be > 0")
    if (!is.numeric(curvature.tol) || curvature.tol <= 0)
        stop("value of 'curvature.tol' must be > 0")
    if (!is.numeric(linesearch.maxit) || linesearch.maxit <= 0)
        stop("maximum number of linesearch iterations must be > 0")

    list(epsilon = epsilon, maxit = maxit, qr.tol = qr.tol,
         improve.tol = improve.tol, curvature.tol = curvature.tol,
         linesearch.method = linesearch.method,
         linesearch.maxit = linesearch.maxit, trace = trace)
}


firthglm.eval <- function(coefficients, x, y, weights, offset, family, control)
{
    ret <- list(indomain=FALSE) # default return value

    # model parameters
    eta <- drop(offset + x %*% coefficients)
    if (!is.null(family$valideta) && !family$valideta(eta))
        return(ret)

    mu <- family$linkinv(eta)
    if (!is.null(family$validmu) && !family$validmu(mu))
        return(ret)

    varmu <- family$variance(mu)
    if (!all(varmu > 0 & is.finite(varmu)))
        return(ret)

    mu.eta <- family$mu.eta(eta)
    skewmu <- family$skewness(mu)

    # qr decomposition
    wt <- weights * mu.eta^2 / varmu
    xw <- x * sqrt(wt)
    qr <- qr(xw, LAPACK=TRUE)
    R <- qr.R(qr)
    q <- qr.Q(qr)
    rownames(R) <- colnames(R)

    # deviance, score, residuals, hat
    dev <- sum(family$dev.resids(y, mu, weights))

    residuals <- (y - mu)/mu.eta
    score <- drop(t(xw) %*% (sqrt(wt) * residuals))
    #hat <- colSums(backsolve(R, t(xw[,qr$pivot]), transpose=TRUE)^2)
    hat <- rowSums(q^2)
    names(hat) <- names(y)

    # penalty
    logdet <- 2 * sum(log(abs(diag(R))))
    logdet.residuals <- ifelse(wt == 0, 0,
                               hat * skewmu * sqrt(varmu)
                               / (weights * mu.eta))
    logdet.grad <- drop(t(xw) %*% (sqrt(wt) * logdet.residuals))

    # modified quantities
    dev.modified <- dev - logdet
    score.modified <- score + 0.5 * logdet.grad
    residuals.modified <- residuals + 0.5 * logdet.residuals

    # domain check
    indomain <- (is.finite(dev.modified)
                 && all(is.finite(score.modified)))

    ## Hessian, sandwich covariance estimator
    ## kurtmu <- family$kurtosis(mu)
    ## kq <- t(q) %*% (ifelse(weights == 0, 0, kurtmu * hat / weights) * q)
    ## q2 <- do.call(cbind, lapply(seq_len(ncol(q)), function(i) q[,i] * q))
    ## sq2 <- ifelse(weights == 0, 0, skewmu / sqrt(weights)) * q2
    ## H <- diag(ncol(q)) - 0.5 * (kq - tcrossprod(t(q) %*% sq2))
    ##
    ## x.modified <- q %*% (H %*% R[,order(qr$pivot),drop=FALSE])
    ## rownames(x.modified) <- rownames(x)
    ## colnames(x.modified) <- colnames(x)
    ##
    ## qr.modified <- qr(x.modified, LAPACK=TRUE)
    ## R.modified <- qr.R(qr.modified)
    ## rownames(R.modified) <- colnames(R.modified)

    list(eta = eta, mu = mu,
         residuals = residuals, residuals.modified = residuals.modified,
         R = R, rank = qr$rank, qr = qr,
         weights = wt, prior.weights = weights,
         deviance = dev, deviance.modified = dev.modified,
         score = score, score.modified = score.modified,
         indomain = indomain)
}


firthglm.fit <-
    function(x, y, weights = rep(1, nobs), start = NULL, etastart = NULL,
             mustart = NULL, offset = rep(0, nobs), family = gaussian(),
             control = list(...), intercept = TRUE, singular.ok = TRUE, ...)
{
    # control
    control <- do.call("firthglm.control", control)

    # design matrix, dimensions
    x <- as.matrix(x)
    xnames <- dimnames(x)[[2L]]
    ynames <- if(is.matrix(y)) rownames(y) else names(y)
    nobs <- NROW(y)
    nvar <- ncol(x)

    # weights, offset
    if (is.null(weights))
        weights <- rep.int(1, nobs)
    if (is.null(offset))
        offset <- rep.int(0, nobs)

    # family
    family <- firthglm.family(family)

    # determine valid range of eta values
    etamax <- .Machine$double.xmax
    etamin <- -(etamax)

    # initial parameters
    n <- NULL # this gets overwritten by eval(family$initizlize)
    if (is.null(mustart)) {
        eval(family$initialize) ## calculates mustart, may change y,
                                ##   weights, and n
    } else {
        mukeep <- mustart
        eval(family$initialize)
        mustart <- mukeep
    }

    eta <- offset
    if (nobs > 0L) {
        mu <- family$linkinv(offset)
        mu.eta <- family$mu.eta(eta)
        varmu <- family$variance(mu)
    } else {
        mu <- numeric()
        varmu <- numeric()
        mu.eta <- numeric()
    }
    wt <- weights * mu.eta^2 / varmu

    qr <- qr(x * sqrt(wt), tol=control$qr.tol)
    rank <- qr$rank
    pivot <- qr$pivot

    # compute initial parameters
    if (!is.null(start)) {
        if (length(start) != nvar) {
            stop(gettextf(paste("length of 'start' should equal %d",
                                "and correspond to initial coefs for %s"),
                          nvar, paste(deparse(xnames), collapse=", ")),
                 domain=NA)
        }
    } else {
        if (is.null(etastart)) {
            if (length(mustart) > 0L) {
                etastart <- family$linkfun(mustart)
            } else {
                etastart <- numeric()
            }
        }
        start <- qr.coef(qr, sqrt(weights) * (etastart - offset))
        start[is.na(start)] <- 0
    }

    # check for rank-deficiency
    xorig <- x
    if (rank < nvar) {
        if (!singular.ok)
            stop("singular fit encountered")
        xdrop <- x[,pivot[(rank+1L):nvar],drop=FALSE]
        x <- x[,pivot[seq_len(rank)],drop=FALSE]
        start <- start[pivot[seq_len(rank)]]
    }

    if (rank == 0L) { # empty model
        coefficients <- rep(NA, nvar)
        names(coefficients) <- xnames
        R <- matrix(NA, 0, 0)
        effects <- qr.qty(qr, sqrt(wt) * y)
        if (!is.null(xnames))
            names(effects) <- character(nobs)

        dev <- sum(family$dev.resids(y, mu, weights))
        nulldev <- dev
        aic <- family$aic(y, n, mu, weights, dev)
        residuals <- (y - mu)/mu.eta

        return(list(coefficients = coefficients,
             residuals = residuals,
             fitted.values = mu, effects = effects,
             R = R, qr = qr, rank = rank,
             family = family, linear.predictors = eta,
             deviance = dev, aic = aic, null.deviance = nulldev,
             iter = 0L, eval = 0L, weights = wt, prior.weights = weights,
             df.residual = nobs, df.null = nobs, y = y, converged = TRUE,
             boundary = FALSE))
    }

    objective <- function(coef)
        firthglm.eval(coef, x, y, weights, offset, family, control)

    coef0 <- start
    obj0 <- objective(coef0)

    if (!obj0$indomain)
        stop("cannot find valid starting values: please specify some",
             call. = FALSE)

    conv <- FALSE
    eval <- 1
    ftol <- control$improve.tol
    gtol <- control$curvature.tol

    for (iter in seq_len(control$maxit)) {
        if (control$trace)
            cat("Penalized deviance = ", obj0$deviance.modified,
                " Iterations - ", iter, "\n", sep = "")

        eta0 <- obj0$eta
        val0 <- 0.5 * (obj0$deviance.modified)
        grad0 <- -(obj0$score.modified)
        search <- numeric(rank)
        search[obj0$qr$pivot] <-
            backsolve(obj0$R, backsolve(obj0$R, transpose=TRUE,
                                        obj0$score.modified[obj0$qr$pivot]))
        # This is mathematically equivalent, but (apparently) less stable:
        #   search <- qr.coef(obj0$qr,
        #                     sqrt(obj0$weights) * obj0$residuals.modified)
        deriv0 <- sum(search * grad0)
        search.eta <- drop(x %*% search)

        # Test for convergence
        if ((deriv0)^2 <= 2 * control$epsilon) {
            conv <- TRUE
            break
        }

        # Fall back to gradient descent if Hessian is ill-conditioned
        if (deriv0 >= 0 || .kappa_tri(obj0$R, LINPACK=FALSE) >= 1e8) {
            search <- -grad0
            deriv0 <- sum(search * grad0)
            search.eta <- drop(x %*% search)
        }

        # determine maximum step size to ensure
        #   |eta[i] - eta0[i]| < 10 * (|eta0[i]| + 1)
        #   for all i
        step.max <- 10 * min((abs(eta0) + 1) / abs(search.eta))
        step.max <- min(step.max, .Machine$double.xmax)

        # determine minimum step size to ensure
        #   |eta[i] - eta0[i]| > eps * (|eta0[i]| + 1)
        #   for at least one i
        step.min <- .Machine$double.eps * min((abs(eta0) + 1) / abs(search.eta))
        step.min <- max(step.min, .Machine$double.xmin)

        # determine initial step, shrinking step.max if necessary
        if (step.min <= 1.0 && 1.0 <= step.max) {
            step0 <- 1.0
        } else {
            step0 <- step.min + 0.5 * (step.max - step.min)
        }
        repeat {
            coef <- coef0 + step0 * search
            obj <- objective(coef)
            eval <- eval + 1
            if (obj$indomain)
                break

            step.max <- step0
            if (step0 < 0.01) {
                step0 <- 2^(0.5 * log2(step.min) + 0.5 * log2(step.max))
            } else {
                step0 <- step.min + 0.5 * (step.max - step.min)
            }
            stopifnot(step0 > step.min)
        }

        # perform line search
        if (control$linesearch.method == "linesearch") {
            lsctrl <- linesearch.control(value.tol = ftol, deriv.tol = gtol,
                                         step.min = step.min,
                                         step.max = step.max)
            ls <- linesearch(val0, deriv0, step0, control = lsctrl)
        } else if (control$linesearch.method == "backtrack") {
            lsctrl <- backtrack.control(value.tol = ftol, step.min = step.min,
                                        contract = 0.5)
            ls <- backtrack(val0, deriv0, step0, control = lsctrl)
        } else if (control$linesearch.method == "blindsearch") {
            lsctrl <- blindsearch.control()
            ls <- blindsearch(val0, deriv0, step0, control = lsctrl)
        } else {
            stop("unrecognized line search method:",
                 control$linesearch.method)
        }

        for (lsiter in seq_len(control$linesearch.maxit)) {
            val <- 0.5 * (obj$deviance.modified)
            grad <- -(obj$score.modified)
            deriv <- sum(search * grad)

            ls <- update(ls, val, deriv)
            if (ls$converged)
                break

            if (control$trace)
                  cat("New step size (", ls$step, ");",
                      " current modified deviance = ",
                      obj$deviance.modified, "\n", sep = "")
            coef <- coef0 + ls$step * search
            obj <- objective(coef)
            eval <- eval + 1
            stopifnot(obj$indomain)
        }

        if (!ls$converged) {
            warning("firthglm.fit: line search failed to converge")
            break
        }

        coef0 <- coef
        obj0 <- obj
        rm(coef, obj)
    }

    eta <- obj0$eta
    mu <- obj0$mu
    dev <- obj0$deviance
    wt <- obj0$weights

    if (rank < nvar) {
        coefficients <- rep(NA, nvar)
        coefficients[pivot[1L:rank]] <- coef0
        names(coefficients) <- xnames
    } else {
        coefficients <- coef0
    }

    # qr.  This is tricky; we can't just call qr(sqrt(wt) * xorig), because
    # we need control of the pivoting
    qr1 <- obj0$qr
    ## qr1 <- obj0$qr.modified
    if (rank < nvar) {
        useLAPACK <- attr(qr1, "useLAPACK")
        i1 <- seq_len(rank)
        i2 <- rank + seq_len(nobs - rank)
        j1 <- seq_len(rank)
        j2 <- (rank+1L):nvar

        x2 <- qr.qty(qr1, sqrt(wt) * xdrop)
        ## q0 <- qr.Q(obj0$qr)
        ## H <- obj0$H
        ## x2 <- qr.qty(qr1, q0 %*% (H %*% (t(q0) %*% (sqrt(wt) * xdrop))))
        x21 <- x2[i1,,drop=FALSE]
        x22 <- x2[i2,,drop=FALSE]

        # LAPACK fails if nrow(x22) == 0
        LAPACK <- !is.null(useLAPACK) && useLAPACK
        if (LAPACK && nrow(x22) == 0L) {
            qr22 <- structure(list(qr = x22, rank = 0L, qraux = numeric(),
                                   pivot = seq_len(ncol(x22))),
                              useLAPACK=TRUE, class="qr")
        } else {
            qr22 <- qr(x22, LAPACK=LAPACK)
        }

        qr <- list()

        qr$qr <- matrix(0, nobs, nvar)
        rownames(qr$qr) <- rownames(qr1$qr)
        colnames(qr$qr) <- c(colnames(qr1$qr), colnames(qr22$qr))
        qr$qr[,j1] <- qr1$qr
        qr$qr[i1,j2] <- x21[,qr22$pivot]
        qr$qr[i2,j2] <- qr22$qr

        qr$rank <- rank
        qr$qraux <- c(qr1$qraux, qr22$qraux)
        qr$pivot <- c(pivot[j1][qr1$pivot], pivot[j2][qr22$pivot])

        class(qr) <- "qr"
        if (!is.null(useLAPACK))
            attr(qr, "useLAPACK") <- useLAPACK ## Important!
    } else {
        qr <- qr1
        qr$pivot <- pivot[qr1$pivot]
    }

    # effects
    effects <- qr.qty(qr, sqrt(wt) * y)
    if (!is.null(xnames))
        names(effects) <- c(xnames[qr$pivot[1L:rank]], rep("", nobs - rank))

    # aic
    aic <- family$aic(y, n, mu, weights, dev) + 2 * rank

    # null deviance
    wtdmu <- if (intercept)
        sum(weights * y)/sum(weights)
    else family$linkinv(offset)
    nulldev <- sum(family$dev.resids(y, wtdmu, weights))

    # df
    n.ok <- nobs - sum(weights == 0)
    nulldf <- n.ok - as.integer(intercept)
    resdf <- n.ok - rank

    list(coefficients = coefficients,
         residuals = obj0$residuals,
         fitted.values = mu, effects = effects,
         R = obj0$R, qr = qr, rank = rank,
         family = family, linear.predictors = eta,
         deviance = dev, aic = aic, null.deviance = nulldev,
         iter = iter, eval = eval, weights = wt, prior.weights = weights,
         df.residual = resdf, df.null = nulldf, y = y, converged = conv,
         boundary = FALSE)
}
