## This file contains all the adjustment functions

#' Adjust for latent factors, after rotationn
#'
#' @param corr.margin marginal correlations, p*d1 matrix
#' @param n sample size
#' @param X.cov estimated second moment of X, d*d matrix
#' @param Gamma estimated confounding effects, p*r matrix
#' @param Sigma diagonal of the estimated noise covariance, p*1 vector
#' @param method adjustment method
#' @param psi derivative of the loss function in robust regression, choices are 
#'			  \code{psi.huber}, \code{psi.bisquare}and \code{psi.hampel}
#' @param nc position of the negative controls
#' @param nc.var.correction correct asymptotic variance based on our formula
#' 
#' @details The function essentially runs a regression of \code{corr.margin} ~ \code{Gamma}. 
#' The sample size \code{n} is needed to have the right scale.
#' 
#' This function should only be called if you know what you are doing. 
#' Most of the time you want to use the main function \code{\link{cate}} to adjust for confounders.
#' 
#' @seealso \code{\link{cate}}
#' 
#' @return a list of objects
#' \describe{
#' \item{alpha}{estimated alpha, r*d1 matrix}
#' \item{beta}{estimated beta, p*d1 matrix}
#' \item{beta.cov.row}{estimated row covariance of \code{beta}, a length p vector}
#' \item{beta.cov.col}{estimated column covariance of \code{beta}, a d1*d1 matrix}
#' }
#'
#' @import stats MASS
#'
#' @export
#'
adjust.latent <- function(corr.margin,
						  n, 
                          X.cov,
                          Gamma,
                          Sigma,
                          method = c("rr", "nc", "lqs"),
                          psi = psi.huber,
                          nc = NULL,
                          nc.var.correction = TRUE) {

    # match arguments
    method <- match.arg(method, c("rr", "nc", "lqs"))

    d <- ncol(X.cov)
    d1 <- ncol(corr.margin)
    p <- nrow(Gamma)
    r <- ncol(Gamma)

    output <- list()
    output$alpha <- matrix(0, r, d1)
    if (method == "rr") {
        for (i in 1:d1) {
            rlm.output <- rlm.cate(x = Gamma, y = corr.margin[, i],
                                  known.scale = sqrt(Sigma)/sqrt(n),
                                   psi = psi)
            output$alpha[, i] <- rlm.output$coefficients
        }
        output$beta <- corr.margin - Gamma %*% output$alpha
        output$beta.cov.row <- Sigma
        output$beta.cov.col <- ginv(X.cov)[(d-d1+1):d, (d-d1+1):d, drop = F] + t(output$alpha) %*% output$alpha
    } else if (method == "nc") { 
        lm.output <- lm(corr.margin[nc, ] ~ Gamma[nc, ] - 1, weights = 1 / Sigma[nc])
        output$alpha <- lm.output$coefficients
        output$alpha <- as.matrix(output$alpha)
        output$beta <- corr.margin - Gamma %*% output$alpha
        if (nc.var.correction) {
            output$beta.cov.row <- Sigma + 
            diag(Gamma %*% solve(t(Gamma[nc,]) %*% diag(1/Sigma[nc]) %*% Gamma[nc, ]) %*% t(Gamma))
        } else {
            output$beta.cov.row <- Sigma
        }
        output$beta.cov.col <- ginv(X.cov)[(d-d1+1):d, (d-d1+1):d, drop = F] + t(output$alpha) %*% output$alpha
    } else if (method == "lqs") {
        for (i in 1:d1) {
            lqs.output <- lqs(Gamma, corr.margin[, i], intercept = FALSE, method = "lts")
            output$alpha[, i] <- lqs.output$coefficients
        }
        output$beta <- corr.margin - Gamma %*% output$alpha
        output$beta.cov.row <- Sigma
        output$beta.cov.col <- ginv(X.cov)[(d-d1+1):d, (d-d1+1):d, drop = F] + t(output$alpha) %*% output$alpha
    }
    colnames(output$alpha) <- colnames(corr.margin)

    return(output)

}

#' Robust linear regression
#'
#' @description This function is slightly modified from rlm.default in the package MASS.
#' The only difference is an option of "known.scale", which is an input vector of the same length of y.
#' 
#' @import stats
#'
#' @keywords internal
#'
rlm.cate <- function (x, y, weights, ..., w = rep(1, nrow(x)), init = "ls",
    psi = psi.huber, scale.est = c("MAD", "Huber", "proposal 2"),
    k2 = 1.345, method = c("M", "MM"), wt.method = c("inv.var",
        "case"), maxit = 100, acc = 1e-04, test.vec = "resid",
    lqs.control = NULL, known.scale = NULL){
    irls.delta <- function(old, new) sqrt(sum((old - new)^2)/max(1e-20,
        sum(old^2)))
    irls.rrxwr <- function(x, w, r) {
        w <- sqrt(w)
        max(abs((matrix(r * w, 1L, length(r)) %*% x)/sqrt(matrix(w,
            1L, length(r)) %*% (x^2))))/sqrt(sum(w * r^2))
    }
    wmad <- function(x, w) {
        o <- sort.list(abs(x))
        x <- abs(x)[o]
        w <- w[o]
        p <- cumsum(w)/sum(w)
        n <- sum(p < 0.5)
        if (p[n + 1L] > 0.5)
            x[n + 1L]/0.6745
        else (x[n + 1L] + x[n + 2L])/(2 * 0.6745)
    }
    method <- match.arg(method)
    wt.method <- match.arg(wt.method)
    nmx <- deparse(substitute(x))
    if (is.null(dim(x))) {
        x <- as.matrix(x)
        colnames(x) <- nmx
    }
    else x <- as.matrix(x)
    if (is.null(colnames(x)))
        colnames(x) <- paste("X", seq(ncol(x)), sep = "")
    if (qr(x)$rank < ncol(x))
        stop("'x' is singular: singular fits are not implemented in 'rlm'")
    if (!(any(test.vec == c("resid", "coef", "w", "NULL")) ||
        is.null(test.vec)))
        stop("invalid 'test.vec'")
    xx <- x
    yy <- y
    if (!missing(weights)) {
        if (length(weights) != nrow(x))
            stop("length of 'weights' must equal number of observations")
        if (any(weights < 0))
            stop("negative 'weights' value")
        if (wt.method == "inv.var") {
            fac <- sqrt(weights)
            y <- y * fac
            x <- x * fac
            wt <- NULL
        }
        else {
            w <- w * weights
            wt <- weights
        }
    }
    else wt <- NULL
    if (method == "M") {
        scale.est <- match.arg(scale.est)
        if (!is.function(psi))
            psi <- get(psi, mode = "function")
        arguments <- list(...)
        if (length(arguments)) {
            pm <- pmatch(names(arguments), names(formals(psi)),
                nomatch = 0L)
            if (any(pm == 0L))
                warning("some of ... do not match")
            pm <- names(arguments)[pm > 0L]
            formals(psi)[pm] <- unlist(arguments[pm])
        }
        if (is.character(init)) {
            temp <- if (init == "ls")
                lm.wfit(x, y, w, method = "qr")
            else if (init == "lts") {
                if (is.null(lqs.control))
                  lqs.control <- list(nsamp = 200L)
                do.call("lqs", c(list(x, y, intercept = FALSE),
                  lqs.control))
            }
            else stop("'init' method is unknown")
            coef <- temp$coefficient
            resid <- temp$residuals
        }
        else {
            if (is.list(init))
                coef <- init$coef
            else coef <- init
            resid <- drop(y - x %*% coef)
        }
    }
    else if (method == "MM") {
        scale.est <- "MM"
        temp <- do.call("lqs", c(list(x, y, intercept = FALSE,
            method = "S", k0 = 1.548), lqs.control))
        coef <- temp$coefficients
        resid <- temp$residuals
        psi <- psi.bisquare
        if (length(arguments <- list(...)))
            if (match("c", names(arguments), nomatch = 0L)) {
                c0 <- arguments$c
                if (c0 > 1.548)
                  formals(psi)$c <- c0
                else warning("'c' must be at least 1.548 and has been ignored")
            }
        scale <- temp$scale
    }
    else stop("'method' is unknown")
    done <- FALSE
    conv <- NULL
    n1 <- (if (is.null(wt))
        nrow(x)
    else sum(wt)) - ncol(x)
    theta <- 2 * pnorm(k2) - 1
    gamma <- theta + k2^2 * (1 - theta) - 2 * k2 * dnorm(k2)
    ###########################################################
    ############Only part changed from rlm.default#############
    ###########################################################
    if (!is.null(known.scale)) {
        scale <- known.scale
    } else {
        if (scale.est != "MM")
            scale <- if (is.null(wt))
                         mad(resid, 0)
                     else wmad(resid, wt)
    }
    ###########################################################
    ############Only part changed from rlm.default#############
    ###########################################################
    for (iiter in 1L:maxit) {
        if (!is.null(test.vec))
            testpv <- get(test.vec)
        if ((scale.est != "MM") && (is.null(known.scale))) {
            ######################################################
            #the second condition aboveis added from rlm.default #
            ######################################################
            scale <- if (scale.est == "MAD")
                if (is.null(wt))
                  median(abs(resid))/0.6745
                else wmad(resid, wt)
            else if (is.null(wt))
                sqrt(sum(pmin(resid^2, (k2 * scale)^2))/(n1 *
                  gamma))
            else sqrt(sum(wt * pmin(resid^2, (k2 * scale)^2))/(n1 *
                gamma))
            if (scale == 0) {
                done <- TRUE
                break
            }
        }
        w <- psi(resid/scale)
        if (!is.null(wt))
            w <- w * weights
        temp <- lm.wfit(x, y, w, method = "qr")
        coef <- temp$coefficients
        resid <- temp$residuals
        if (!is.null(test.vec))
            convi <- irls.delta(testpv, get(test.vec))
        else convi <- irls.rrxwr(x, w, resid)
        conv <- c(conv, convi)
        done <- (convi <= acc)
        if (done)
            break
    }
    if (!done)
        warning(gettextf("'rlm' failed to converge in %d steps",
            maxit), domain = NA)
    fitted <- drop(xx %*% coef)
    cl <- match.call()
    cl[[1L]] <- as.name("rlm")
    fit <- list(coefficients = coef, residuals = yy - fitted,
        wresid = resid, effects = temp$effects, rank = temp$rank,
        fitted.values = fitted, assign = temp$assign, qr = temp$qr,
        df.residual = NA, w = w, s = scale, psi = psi, k2 = k2,
        weights = if (!missing(weights)) weights, conv = conv,
        converged = done, x = xx, call = cl)
    class(fit) <- c("rlm", "lm")
    fit
}

## #' Modified from lqs.default
## lqs.cate <- function (x, y, intercept = TRUE, method = c("lts", "lqs", "lms", "S"), quantile, control = lqs.control(...), k0 = 1.548, seed,
##           ###added known.scale
##           known.scale = NULL, ...)
## {
##     lqs.control <- function(psamp = NA, nsamp = "best", adjust = TRUE) list(psamp = psamp,
##         nsamp = nsamp, adjust = adjust)
##     n <- length(y)
##     nmx <- deparse(substitute(x))
##     if (is.null(dim(x))) {
##         x <- as.matrix(x)
##         colnames(x) <- nmx
##     }
##     else x <- as.matrix(x)
##     p <- ncol(x)
##     if (any(is.na(x)) || any(is.na(y)))
##         stop("missing values are not allowed")
##     nm <- colnames(x)
##     if (is.null(nm))
##         nm <- if (p > 1)
##             paste("X", 1L:p, sep = "")
##         else if (p == 1)
##             "X"
##         else NULL
##     if (intercept) {
##         x <- cbind(1, x)
##         nm <- c("(Intercept)", nm)
##     }
##     p <- ncol(x)
##     if (nrow(x) != n)
##         stop("'x' and 'y' must have the same number of rows")
##     method <- match.arg(method)
##     lts <- 0
##     beta <- 0
##     if (method == "lqs" && missing(quantile))
##         quantile <- floor((n + p + 1)/2)
##     if (method == "lms")
##         quantile <- floor((n + 1)/2)
##     if (method == "lts") {
##         lts <- 1
##         if (missing(quantile))
##             quantile <- floor(n/2) + floor((p + 1)/2)
##     }
##     if (method == "S") {
##         lts <- 2
##         beta <- 0.5
##         quantile <- ceiling(n/2)
##         chi <- function(u, k0) {
##             u <- (u/k0)^2
##             ifelse(u < 1, 3 * u - 3 * u^2 + u^3, 1)
##         }
##     }
##     if (quantile > n - 1)
##         stop(gettextf("'quantile' must be at most %d", n - 1),
##             domain = NA)
##     ps <- control$psamp
##     if (is.na(ps))
##         ps <- p
##     if (ps < p) {
##         ps <- p
##         warning("'ps' must be at least 'p'")
##     }
##     adj <- control$adjust & intercept
##     nsamp <- eval(control$nsamp)
##     nexact <- choose(n, ps)
##     if (is.character(nsamp) && nsamp == "best") {
##         nsamp <- if (nexact < 5000)
##             "exact"
##         else "sample"
##     }
##     else if (is.numeric(nsamp) && nsamp > nexact) {
##         warning(sprintf(ngettext(nexact, "only %d set, so all sets will be tried",
##             "only %d sets, so all sets will be tried"), nexact),
##             domain = NA)
##         nsamp <- "exact"
##     }
##     samp <- nsamp != "exact"
##     if (samp) {
##         if (nsamp == "sample")
##             nsamp <- min(500 * ps, 3000)
##     }
##     else nsamp <- nexact
##     if (samp && !missing(seed)) {
##         if (exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) {
##             seed.keep <- get(".Random.seed", envir = .GlobalEnv,
##                 inherits = FALSE)
##             on.exit(assign(".Random.seed", seed.keep, envir = .GlobalEnv))
##         }
##         assign(".Random.seed", seed, envir = .GlobalEnv)
##     }
##     z <- .C(lqs_fitlots, as.double(x), as.double(y), as.integer(n),
##         as.integer(p), as.integer(quantile), as.integer(lts),
##         as.integer(adj), as.integer(samp), as.integer(ps), as.integer(nsamp),
##         crit = double(1), sing = integer(1L), bestone = integer(ps),
##         coefficients = double(p), as.double(k0), as.double(beta))[c("crit",
##         "sing", "coefficients", "bestone")]
##     if (z$sing == nsamp)
##         stop("'lqs' failed: all the samples were singular", call. = FALSE)
##     z$sing <- paste(z$sing, "singular samples of size", ps, "out of",
##         nsamp)
##     z$bestone <- sort(z$bestone)
##     names(z$coefficients) <- nm
##     fitted <- drop(x %*% z$coefficients)
##     z$fitted.values <- fitted
##     z$residuals <- y - fitted
##     c1 <- 1/qnorm((n + quantile)/(2 * n))
##     s <- if (lts == 1)
##         sqrt(z$crit/quantile)/sqrt(1 - 2 * n * dnorm(1/c1)/(quantile *
##             c1))
##     else if (lts == 0)
##         sqrt(z$crit) * c1
##     else z$crit
##     res <- z$residuals
##     ind <- abs(res) <= 2.5 * s
##     s2 <- sum(res[ind]^2)/(sum(ind) - p)
##     z$scale <- c(s, sqrt(s2))
##     if (method == "S") {
##         psi <- function(u, k0) (1 - pmin(1, abs(u/k0))^2)^2
##         resid <- z$residuals
##         scale <- s
##         for (i in 1L:30) {
##             w <- psi(resid/scale, k0)
##             temp <- lm.wfit(x, y, w, method = "qr")
##             resid <- temp$residuals
##             s2 <- scale * sqrt(sum(chi(resid/scale, k0))/((n -
##                 p) * beta))
##             if (abs(s2/scale - 1) < 1e-05)
##                 break
##             scale <- s2
##         }
##         z$coef <- temp$coefficients
##         z$fitted.values <- temp$fitted.values
##         z$residuals <- resid
##         z$scale <- scale
##     }
##     class(z) <- "lqs"
##     z
## }
