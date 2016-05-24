
### <FIXME> rename fixed to coef and allow for specification of coefs, 
###         ie fitted models? </FIXME>
.mlt_setup <- function(model, data, y, offset = NULL, fixed = NULL) {

    response <- model$response
    stopifnot(length(response) == 1)
    todistr <- model$todistr

    eY <- .mm_exact(model, data = data, response = response, object = y)
    iY <- .mm_interval(model, data = data, response = response, object = y)

    if (is.null(eY)) {
        Y <- iY$Yleft
    } else {
        Y <- eY$Y
    }

    ui <- attr(Y, "constraint")$ui
    ci <- attr(Y, "constraint")$ci

    if (!is.null(fixed)) {
        stopifnot(all(names(fixed) %in% colnames(Y)))
        fix <- colnames(Y) %in% names(fixed)
        ui <- ui[,!fix,drop = FALSE]
        .parm <- function(beta) {
            ret <- numeric(ncol(Y))
            ret[fix] <- fixed
            ret[!fix] <- beta
            ret
        }
    } else {
        .parm <- function(beta) beta
        fix <- rep(FALSE, ncol(Y))
    } 

    exact <- .exact(y)

    ll <- function(beta) {
        ret <- numeric(nrow(data))
        if (any(exact))
            ret[exact] <- .mlt_loglik_exact(todistr, eY$Y, eY$Yprime, offset[exact], eY$trunc)(.parm(beta))
        ret[!exact] <- .mlt_loglik_interval(todistr, iY$Yleft, iY$Yright, offset[!exact], iY$trunc)(.parm(beta))
        return(ret)
    }
    sc <- function(beta) {
        ret <- matrix(0, nrow = nrow(data), ncol = length(fix))
        if (any(exact))
            ret[exact,] <- .mlt_score_exact(todistr, eY$Y, eY$Yprime, offset[exact], eY$trunc)(.parm(beta))
        ret[!exact,] <- .mlt_score_interval(todistr, iY$Yleft, iY$Yright, offset[!exact], iY$trunc)(.parm(beta))
        return(ret[, !fix, drop = FALSE])
    }
    he <- function(beta, weights) {
        ret <- 0
        if (any(exact))
            ret <- ret + .mlt_hessian_exact(todistr, eY$Y, eY$Yprime, offset[exact], eY$trunc, weights[exact])(.parm(beta))
        if (any(!exact))
            ret <- ret + .mlt_hessian_interval(todistr, iY$Yleft, iY$Yright, offset[!exact], iY$trunc, weights[!exact])(.parm(beta))
        return(ret[!fix, !fix, drop = FALSE])
    }

    loglikfct <- function(beta, weights) -sum(weights * ll(beta))
    score <- function(beta, weights) weights * sc(beta)
    scorefct <- function(beta, weights) -colSums(score(beta, weights))

    if (all(!is.finite(ci))) {
        ui <- ci <- NULL
    } else {
        ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
        ci <- ci[is.finite(ci)]
        r0 <- rowSums(abs(ui)) == 0
        ui <- ui[!r0,,drop = FALSE]
        ci <- ci[!r0]
        if (nrow(ui) == 0) ui <- ci <- NULL
        ### ci <- ci + sqrt(.Machine$double.eps) ### we need ui %*% theta > ci, not >= ci
    }

    optimfct <- function(theta, weights, scale = FALSE, quiet = TRUE, ...) {
        control <- list(...)
        if (scale) {
            Ytmp <- Y
            Ytmp[!is.finite(Ytmp)] <- NA
            sc <- apply(abs(Ytmp[, !fix, drop = FALSE]), 2, max, na.rm = TRUE)
            lt1 <- sc < 1.1
            gt1 <- sc >= 1.1
            sc[gt1] <- 1 / sc[gt1]
            sc[lt1] <- 1
            f <- function(gamma) loglikfct(sc * gamma, weights)
            g <- function(gamma) scorefct(sc * gamma, weights) * sc
            theta <- theta / sc
            if (!is.null(ui))
                ui <- t(t(ui) * sc)
        } else {
            f <- function(gamma) loglikfct(gamma, weights)
            g <- function(gamma) scorefct(gamma, weights)
        }
        if (!is.null(ui)) {
            ret <- BBoptim(par = theta, fn = f, gr = g, project = "projectLinear",
                           projectArgs = list(A = ui, b = ci, meq = 0), control = control, 
                           quiet = quiet)
        } else {
            ret <- BBoptim(par = theta, fn = f, gr = g, control = control, quiet = quiet)
        }
        if (quiet & (ret$convergence != 0))
            warning("Optimisation did not converge")
        ### degrees of freedom: number of free parameters, ie #parm NOT meeting the constraints
        ret$df <- length(ret$par)
        ### <FIXME> check on alternative degrees of freedom
#        if (!is.null(ui)) 
#            ret$df <- ret$df - sum(ui %*% ret$par - ci < .Machine$double.eps)
        ### </FIXME>
        if (scale) ret$par <- ret$par * sc
        return(ret)
    }

    coef <- rep(NA, length(fix))
    coef[fix] <- fixed
    names(coef) <- colnames(Y)

    ret <- list()
    ret$parm <- .parm
    ret$coef <- coef
    ret$fixed <- fixed
    ret$model <- model
    ret$data <- data
    ret$offset <- offset
    ret$response <- response
    ret$todistr <- todistr
    ret$loglik <- loglikfct
    ret$score <- score
    ret$hessian <- he
    ret$optimfct <- optimfct
    class(ret) <- c("mlt_setup", "mlt")
    return(ret)
}

.mlt_start <- function(model, data, y, pstart, offset = NULL, fixed = NULL) {

    stopifnot(length(pstart) == nrow(data))
    if (is.null(offset)) offset <- rep(0, nrow(data))

    response <- model$response
    stopifnot(length(response) == 1)

    eY <- .mm_exact(model, data = data, response = response, object = y)
    iY <- .mm_interval(model, data = data, response = response, object = y)

    if (is.null(eY)) {
        Y <- iY$Yleft
    } else {
        Y <- eY$Y
    }

    ui <- as.matrix(attr(Y, "constraint")$ui)
    ci <- attr(Y, "constraint")$ci

    if (!is.null(fixed)) {
        stopifnot(all(names(fixed) %in% colnames(Y)))
        fix <- colnames(Y) %in% names(fixed)
        ui <- ui[,!fix,drop = FALSE]
    } else {
        fix <- rep(FALSE, ncol(Y))
    } 

    X <- matrix(0, nrow = NROW(y), ncol = ncol(Y))
    if (!is.null(eY))
        X[eY$which,] <- eY$Y
    if (!is.null(iY))
        X[iY$which,] <- iY$Yright
    X[!is.finite(X[,1]),] <- 0

    if (any(fix)) {
        offset <- X[, fix, drop = FALSE] %*% fixed
        X <- X[, !fix, drop = FALSE]
    }

    todistr <- model$todistr
    Z <- todistr$q(pmax(.01, pmin(pstart, .99))) - offset

    dvec <- crossprod(X, Z)
    Dmat <- crossprod(X)
    diag(Dmat) <- diag(Dmat) + 1e-08

    if (all(!is.finite(ci))) {
        ui <- ci <- NULL
    } else {
        ui <- as(ui[is.finite(ci),,drop = FALSE], "matrix")
        ci <- ci[is.finite(ci)]
        r0 <- rowSums(abs(ui)) == 0
        ui <- ui[!r0,,drop = FALSE]
        ci <- ci[!r0]
        if (nrow(ui) == 0) ui <- ci <- NULL
        ci <- ci + sqrt(.Machine$double.eps) ### we need ui %*% theta > ci, not >= ci
    }

    if (!is.null(ui)) {    
        ret <- solve.QP(Dmat, dvec, t(ui), ci, meq = 0)$solution
    } else {
        ret <- lm.fit(x = X, y = Z)$coef
    }
    ret
}

.mlt_fit <- function(object, weights, theta = NULL, scale = FALSE, trace = FALSE, 
                     quiet = TRUE, ...) {

    if (is.null(theta))
        stop(sQuote("mlt"), "needs suitable starting values")

    ### BBoptim issues a warning in case of unsuccessful convergence
    ret <- try(object$optimfct(theta, weights = weights, 
                               trace = trace, scale = scale, quiet = quiet, ...))    

    cls <- class(object)
    object[names(ret)] <- NULL
    object <- c(object, ret)
    object$coef[] <- object$parm(ret$par) ### [] preserves names
    object$theta <- theta ### starting value
    object$scale <- scale ### scaling yes/no
    object$weights <- weights
    object$trace <- trace
    object$quiet <- quiet
    class(object) <- c("mlt_fit", cls)
    
    return(object)
}

mlt <- function(model, data, weights = NULL, offset = NULL, fixed = NULL,
                theta = NULL, pstart = NULL, scale = FALSE,
                checkGrad = FALSE, trace = FALSE, quiet = TRUE, dofit = TRUE, ...) {

    vars <- as.vars(model)
    response <- model$response
    responsevar <- vars[[response]]
    bounds <- bounds(responsevar)
    stopifnot(length(response) == 1)
    y <- R(object = data[[response]])

    if (is.null(weights)) weights <- rep(1, nrow(data))
    if (is.null(offset)) offset <- rep(0, nrow(data))
    stopifnot(nrow(data) == length(weights))
    stopifnot(nrow(data) == length(offset))

    s <- .mlt_setup(model = model, data = data, y = y, 
                    offset = offset, fixed = fixed) 
    if (!dofit) return(s)

    if (is.null(theta)) {
        ### unconditional ECDF, essentially
###        if (is.null(pstart)) pstart <- y$rank / max(y$rank)
        if (is.null(pstart)) pstart <- attr(y, "prob")(weights)(y$approxy) ### y$rank / max(y$rank)
        theta <- .mlt_start(model = model, data = data, y = y, 
                            pstart = pstart, offset = offset, fixed = fixed)
    }

    args <- list(...)
    args$object <- s
    args$weights <- weights
    args$theta <- theta
    args$scale <- scale
    args$trace <- trace
    args$checkGrad <- checkGrad
    args$quiet <- quiet
    ret <- do.call(".mlt_fit", args)
    ret$call <- match.call()
    ret$bounds <- bounds
    ret
}

update.mlt_fit <- function(object, weights, theta, ...) {

    stopifnot(length(weights) == NROW(object$data))
    args <- list(...)
    if (inherits(object, "mlt_fit")) 
        class(object) <- class(object)[-1L]
    args$object <- object
    if (missing(weights)) {
        args$weights <- weights(object)
    } else {
        args$weights <- weights
    }
    if (missing(theta)) {
        args$theta <- object$theta
    } else {
        args$theta <- theta
    }
    args$scale <- object$scale
    args$trace <- object$trace
    args$checkGrad <- object$checkGrad
    args$quiet <- object$quiet
    ret <- do.call(".mlt_fit", args)
    ret$call <- match.call()
    ret
}
