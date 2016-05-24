glmmPQL.fit <- function(X, y, Z,  weights = rep(1, NROW(y)), start = NULL,
                        etastart = NULL, mustart = NULL,
                        offset = rep(0, NROW(y)), family = gaussian(),
                        control = glmmPQL.control(...),
                        sigma = NULL, sigma.fixed = FALSE, ...) {
    matchCall <- as.list(match.call(expand.dots = FALSE))
    dots <- names(matchCall[["..."]])
    dots <- intersect(dots, names(formals(glm)))
    fit0 <- do.call("glm.fit", c(list(X, y, weights, start = start,
                                      etastart = etastart, mustart = mustart,
                                      offset = offset, family = family),
                                 matchCall[dots]))
    w <- fit0$prior.weights
    # QR missing from glm.fit if ncol(X) = 0
    QR <- qr(X)
    R <- qr.R(QR)
    rank <- QR$rank
    p <- ncol(R)
    nm <- colnames(R)[seq(length = rank)]
    if (rank < p) {
        X0 <- X[,colnames(R)[-seq(length = rank)]]
        X <- X[, nm]
    }
    empty <- !length(X)
    if (empty) {
        alpha <- numeric(0)
        Xa <- matrix(0, length(y), 1)
    }
    eta <- fit0$linear.predictors
    residuals <- fit0$residuals
    Y <- eta + residuals - offset #working response
    wy <- fit0$weights # iterative weights
    wY <- sqrt(wy) * Y
    wZ <- sqrt(wy) * Z
    ZWy <- crossprod(wZ, wY)
    ZWZ <- crossprod(wZ, wZ)
    if (!empty) {
        wX <- sqrt(wy) * X
        XWy <- crossprod(wX, wY)
        XWX <- crossprod(wX, wX)
        ZWX <- crossprod(wZ, wX)
        E <- chol(XWX)
        F <- backsolve(E, t(ZWX), transpose = TRUE)
        f <- backsolve(E, XWy, transpose = TRUE)
        ZSy <- ZWy - crossprod(F, f)
        ZSZ <- ZWZ - crossprod(F, F)
    }
    if (is.null(sigma)) sigma <- 0.1
    logtheta <- log(sigma^2)
    conv <- FALSE
    for (i in 1:control$maxiter) {
        ## Update coefficients
        for (j in 1:control$IWLSiter) {
            IZWZD <- ZWZ * sigma^2
            diag(IZWZD) <- 1 + diag(IZWZD)
            A <- chol(IZWZD)
            if (!empty) {
                IZSZD <- ZSZ * sigma^2
                diag(IZSZD) <- 1 + diag(IZSZD)
                G <- chol(IZSZD)
                g <- backsolve(G, ZSy, transpose = TRUE)
                v <- backsolve(G, g)

                B <- backsolve(A, sigma * ZWX, transpose = TRUE)
                C <- chol(XWX - crossprod(B, B))
                b <- backsolve(A, sigma * ZWy, transpose = TRUE)
                c <- backsolve(C, XWy - t(B) %*% b, transpose = TRUE)

                alpha <- backsolve(C, c)
                Xa <- X %*% alpha
                beta <- sigma^2 * v
            }
            else {
                g <- backsolve(A, ZWy, transpose = TRUE)
                v <- backsolve(A, g)

                beta <- sigma^2 * v
            }

            eta <- c(Xa + Z %*% beta + offset)

            ## Update working response & weights
            mu <- family$linkinv(eta)
            mu.eta.val <- family$mu.eta(eta)
            residuals <- (fit0$y - mu)/mu.eta.val
            Y <- eta + residuals - offset
            wy <- w * mu.eta.val^2/family$variance(mu)
            wY <- sqrt(wy) * Y
            wZ <- sqrt(wy) * Z
            ZWy <- crossprod(wZ, wY)
            ZWZ <- crossprod(wZ, wZ)

            if (!empty) {
                wX <- sqrt(wy) * X
                XWy <- crossprod(wX, wY)
                XWX <- crossprod(wX, wX)
                ZWX <- crossprod(wZ, wX)
                E <- chol(XWX)
                F <- backsolve(E, t(ZWX), transpose = TRUE)
                f <- backsolve(E, XWy, transpose = TRUE)
                ZSy <- ZWy - crossprod(F, f)
                ZSZ <- ZWZ - crossprod(F, F)

                score <- c(crossprod(X, wy * residuals),
                           crossprod(Z, wy * residuals) - v)
                diagInfo <- c(diag(XWX), diag(ZWZ))

                if (all(diagInfo < 1e-20) ||
                    all(abs(score) < control$tol * sqrt(control$tol + diagInfo))) {
                    if (sigma.fixed) conv <- TRUE
                    break
                }
            }
            else {
                score <- crossprod(Z, wy * residuals) - v
                diagInfo <- diag(ZWZ)
                if (all(diagInfo < 1e-20) ||
                    all(abs(score) < control$tol * sqrt(control$tol + diagInfo))) {
                    if (sigma.fixed) conv <- TRUE
                    break
                }
            }
        }

        if (!sigma.fixed){
            ## Update sigma
            ## sigma^2 = exp(logtheta)
            ## One Fisher scoring iteration
            IZWZD <- ZWZ * sigma^2
            diag(IZWZD) <- 1 + diag(IZWZD)
            A <- chol(IZWZD)
            if (!empty) {
                IZSZD <- ZSZ * sigma^2
                diag(IZSZD) <- 1 + diag(IZSZD)
                G <- chol(IZSZD)
                g <- backsolve(G, ZSy, transpose = TRUE)
                v <- backsolve(G, g)
                h <- backsolve(G, ZSZ, transpose = TRUE)
                H <- backsolve(G, h)
            }
            else {
                g <- backsolve(A, ZWy, transpose = TRUE)
                v <- backsolve(A, g)
                h <- backsolve(A, ZWZ, transpose = TRUE)
                H <- backsolve(A, h)
            }

            ## Harville p326
            score <- drop(-0.5 * sum(diag(H)) + 0.5 * crossprod(v, v)) *
                sigma^2
            Info <- 0.5 * sum(H^2) * sigma^4

             if (control$trace) {
                 ##B & C eq 5 - still not consistently increasing
                 cat("Iteration ", i,
                     ". Score = ", abs(score) ,
                     "\n", sep = "")
                 flush.console()
            }

            ## check for overall convergence
            if (Info < 1e-20 ||
                abs(score) < control$tol * sqrt(control$tol + Info)){
                conv <- TRUE
                break
            }

            ## Cannot use beta to update t(YXa) %*% Vinv %*% YXa
            ZWYXa <- crossprod(wZ, sqrt(wy) * (Y - Xa))
            optfun <- function(logtheta) {
                IZWZD <- ZWZ * exp(logtheta)
                diag(IZWZD) <- 1 + diag(IZWZD)
                A <- chol(IZWZD)
                if (!empty) {
                    IZSZD <- ZSZ * exp(logtheta)
                    diag(IZSZD) <- 1 + diag(IZSZD)
                    G <- chol(IZSZD)
                    d <- backsolve(A, sqrt(exp(logtheta)) * ZWYXa, transpose = TRUE)
                    sum(log(diag(G))) - 0.5 * crossprod(d, d)
                }
                else {
                    d <- backsolve(A, sqrt(exp(logtheta)) * ZWy, transpose = TRUE)
                    sum(log(diag(A))) - 0.5 * crossprod(d, d)
                }
            }
            optres <- optimize(optfun, c(-10, 10))
            if (optfun(-10) < optfun(optres$minimum))
                sigma <- 0
            else {
                if (abs(optres$minimum - (logtheta + score/Info)) > 0.1)
                    logtheta <- optres$minimum
                else
                    logtheta <- logtheta + score/Info
                sigma <- sqrt(exp(logtheta))
            }
        }
        else if (conv)
            break
    }
    if (!empty) varFix <- chol2inv(C)
    else varFix <- matrix(, 0, 0)
    rownames(varFix) <- colnames(varFix) <- colnames(X)
    fit0$coef[nm] <- alpha
    if (!sigma.fixed)
        varSigma <- sigma^2/(4 * Info)
    else
        varSigma <- NA
    glm <- identical(sigma, 0)
    if (!empty) {
        if (rank < p) QR <- qr(cbind(wX, sqrt(w) * X0))
        else QR <- qr(wX)
        R <- qr.R(QR)
    }
    list(coefficients = structure(fit0$coef, random = beta),
         residuals = residuals,
         fitted.values = mu,
         #effect = ?
         R = if (!empty) R,
         rank = rank,
         qr = if (!empty) QR,
         family = family,
         linear.predictors = eta,
         deviance = if (glm) sum(family$dev.resids(y, mu, w)),
         aic = if (glm)
         family$aic(y, length(y), mu, w, sum(family$dev.resids(y, mu, w))) + 2 * rank,
         null.deviance = if (glm) {
             wtdmu <- family$linkinv(offset)
             sum(family$dev.resids(y, wtdmu, w))
         },
         iter = ifelse(glm, NA, i),
         weights = wy,
         prior.weights = w,
         df.residual = length(y) - rank,
         df.null = if (glm) length(y) - sum(w == 0),
         y = y,
         sigma = sigma, sigma.fixed = sigma.fixed,
         varFix = varFix, varSigma = varSigma, converged = conv)
}
