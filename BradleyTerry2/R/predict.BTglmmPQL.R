predict.BTglmmPQL <- function (object, newdata = NULL, newrandom = NULL,
                             level = 1, type = c("link", "response", "terms"),
                             se.fit = FALSE, terms = NULL,
                             na.action = na.pass, ...) {
    ## only pass on if a glm
    if (object$sigma == 0) {
        if (level != 0) warning("Fixed effects model: setting level to 0")
        return(NextMethod())
    }
    if (!all(level %in% c(0, 1))) stop("Only level %in% c(0, 1) allowed")
    type <- match.arg(type)
    if (!is.null(newdata) || type == "terms") tt <- terms(object)
    if (!is.null(newdata)) {
        ## newdata should give variables in terms formula
        Terms <- delete.response(tt)
        m <- model.frame(Terms, newdata, na.action = na.action,
                         xlev = object$xlevels)
        na.action <- attr(m, "na.action")
        if (!is.null(cl <- attr(Terms, "dataClasses")))
            .checkMFClasses(cl, m)
        D <- model.matrix(Terms, m, contrasts.arg = object$contrasts)
        np <- nrow(D) # n predictions
        offset <- rep(0, np)
        if (!is.null(off.num <- attr(tt, "offset")))
            for (i in off.num) offset <- offset + eval(attr(tt,
                "variables")[[i + 1]], newdata)
        if (!is.null(object$call$offset))
            offset <- offset + eval(object$call$offset, newdata)
    }
    else {
        D <- model.matrix(object)
        newrandom <- object$random
        na.action <- object$na.action
        offset <- object$offset
    }
    if (se.fit == TRUE) {
        sigma <- object$sigma
        w <- sqrt(object$weights)
        if (!is.null(newdata)) {
            wX <- w * model.matrix(object)
            wZ <- w * object$random
        }
        else {
            wX <- w * D
            wZ <- w * newrandom
        }
        XWX <- crossprod(wX)
        XWZ <- crossprod(wX, wZ)
        ZWZ <- crossprod(wZ, wZ)
        diag(ZWZ) <- diag(ZWZ) + 1/sigma^2
        C <- cbind(XWX, XWZ)
        C <- chol(rbind(C, cbind(t(XWZ), ZWZ)))
        if (type == "terms" || level == 0){
            ## work out (chol of inverse of) topleft of C-inv directly
            A <- backsolve(chol(ZWZ), t(XWZ), transpose = TRUE)
            A <- chol(XWX - t(A) %*% A)
        }
    }
    if (type == "terms") { # ignore level
        if (1 %in% level)
            warning("type = \"terms\": setting level to 0", call. = FALSE)
        coef <- coef(object) #fixef
        aa <- attr(D, "assign")
        ll <- attr(tt, "term.labels")
        if (!is.null(terms)) {
            include <- ll %in% terms
            ll <- ll[include]
        }
        hasintercept <- attr(tt, "intercept") > 0L
        if (hasintercept) {
            avx <- colMeans(model.matrix(object))
            termsconst <- sum(avx * coef) #NA coefs?
            D <- sweep(D, 2, avx)
        }
        pred0 <- matrix(ncol = length(ll), nrow = NROW(D))
        colnames(pred0) <- ll
        if (se.fit) {
            A <- chol2inv(A)
            se.pred0 <- pred0
        }
        for (i in seq(length.out = length(ll))){
            ind <- aa == which(attr(tt, "term.labels") == ll[i])
            pred0[, i] <- D[, ind, drop = FALSE] %*% coef[ind]
            if (se.fit) {
                se.pred0[, i] <-  sqrt(diag(D[, ind] %*%
                                            tcrossprod(A[ind, ind], D[, ind])))
            }
        }
        if (hasintercept) attr(pred0, "constant") <- termsconst
        if (se.fit) return(list(fit = pred0, se.fit = se.pred0))
        return(pred0)
    }
    if (0 %in% level) {
        pred0 <- napredict(na.action, c(D %*% coef(object)) + offset)
        if (type == "response")
            pred0 <- family(object)$linkinv(pred0)
        if (se.fit == TRUE) {
            na.act <- attr(na.exclude(pred0), "na.action")
            H <- backsolve(A, t(na.exclude(D)), transpose = TRUE)
            ## se.pred0 <- sqrt(diag(D %*% chol2inv(C)[1:ncol(D), 1:ncol(D)] %*% t(D)))
            se.pred0 <- napredict(na.action, napredict(na.act, sqrt(colSums(H^2))))
           if (type == "response")
               se.pred0 <- se.pred0*abs(family(object)$mu.eta(pred0))
            pred0 <- list(fit = pred0, se.fit = se.pred0)
        }
        if (identical(level, 0)) return(pred0)
    }

    r <- nrow(D)
    ## newrandom should give new design matrix for original random effects
    if (!is.null(newdata)){
        if(is.null(newrandom))
            stop("newdata specified without newrandom")
        if (!is.null(na.action))
            newrandom <- newrandom[-na.action, , drop = FALSE]
    }
    if (!identical(dim(newrandom), c(r, ncol(object$random))))
        stop("newrandom should have ", r, " rows and ",
             ncol(object$random), " columns")
    D <- cbind(D, newrandom)
    coef <- c(coef(object), attr(coef(object), "random"))
    pred <- napredict(na.action, c(D %*% coef) + offset)
    if (type == "response")
        pred <- family(object)$linkinv(pred)
    if (se.fit == TRUE) {
        ##se.pred <- sqrt(diag(D %*% chol2inv(C) %*% t(D)))
        na.act <- attr(na.exclude(pred), "na.action")
        H <- backsolve(C, t(na.exclude(D)), transpose = TRUE)
        se.pred <- napredict(na.action, napredict(na.act, sqrt(colSums(H^2))))
        if (type == "response")
            se.pred <- se.pred*abs(family(object)$mu.eta(pred))
        pred <- list(fit = pred, se.fit = se.pred)
    }

    if (0 %in% level)
        list(population = pred0, individual = pred)
    else pred
}
