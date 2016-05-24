### widekernelpls.fit.R: Kernel PLS fit algorithm for wide data.
### $Id: widekernelpls.fit.R 197 2011-11-15 08:55:40Z bhm $
###
### Implements an adapted version of the algorithm described in
###  Rannar, S., Lindgren, F., Geladi, P. and Wold, S. (1994) A PLS
###  Kernel Algorithm for Data Sets with Many Variables and Fewer
###  Objects.  Part 1: Theory and Algorithm.
###  \emph{Journal of Chemometrics}, \bold{8}, 111--125.

widekernelpls_NA.fit <- function(X, Y, ncomp, stripped = FALSE,
                                 tol = .Machine$double.eps^0.5,
                                 maxit = 100, EM.tol = .Machine$double.eps^0.5,
                                 EM.maxit = 10,
                                 EM.init = c("col", "meancolrow", "colrow",
                                 "colrowmean", "neigh1", "neigh2"),
                                 ...)
{
    ## Initialise
    Y <- as.matrix(Y)
    if(!stripped) {
        ## Save dimnames:
        dnX <- dimnames(X)
        dnY <- dimnames(Y)
    }
    ## Remove dimnames during calculation.
    dimnames(X) <- dimnames(Y) <- NULL

    nobj <- dim(X)[1]
    npred <- dim(X)[2]
    nresp <- dim(Y)[2]

    TT <- U <- matrix(0, ncol = ncomp, nrow = nobj)# scores
    B <- array(0, c(npred, nresp, ncomp))
    In <- diag(nobj)
    nits <- numeric(ncomp)              # for debugging
    if (!stripped) {
        fitted <- array(0, dim = c(nobj, nresp, ncomp))
        Xresvar <- numeric(ncomp)
    }

    ## EM Step 1 (Initialisation)
    EM.init <- match.arg(EM.init)
    Xmissing <- which(is.na(X), arr.ind = TRUE)
    XRmeans <- rowMeans(X, na.rm = TRUE)
    XCmeans <- colMeans(X, na.rm = TRUE)
    Xmean <- mean(X, na.rm = TRUE)
    for (ind in seq_len(nrow(Xmissing))) {
        row <- Xmissing[ind,1]
        col <- Xmissing[ind,2]
        X[row,col] <-
            switch(EM.init,
                   col = XCmeans[col],
                   meancolrow = (XRmeans[row] + XCmeans[col])/2,
                   colrow = XRmeans[row] + XCmeans[col],
                   colrowmean = XRmeans[row] + XCmeans[col] - Xmean,
                   neigh1 = mean(X[row,][col + (-1):1], na.rm = TRUE),
                   neigh2 = mean(X[row,][col + (-2):2], na.rm = TRUE)
                   )
    }
    ## FIXME: make it a function, then do Y as well

    ## EM: Repeat updating estimates
    P.old <- matrix(0, nrow = npred, ncol = ncomp) # Initial values
    EM.nit <- 0                                    # for debugging
    changeX <- change <- numeric(EM.maxit)
    Xs <- Ps <- list()
    repeat {
        EM.nit <- EM.nit + 1            # for debugging
        Xold <- X[Xmissing]    # Save the old estimates for convergence checks

        ## Center variables:
        Xmeans <- colMeans(X)
        X0 <- X - rep(Xmeans, each = nobj)
        Ymeans <- colMeans(Y)
        Y0 <- Y - rep(Ymeans, each = nobj)

        XXt <- tcrossprod(X0)
        YYt <- tcrossprod(Y0)

        if (!stripped) Xtotvar <- sum(diag(XXt))

        for (a in 1:ncomp) {
            XXtYYt <- XXt %*% YYt
            ## This avoids problems with negative eigenvalues due to roundoff
            ## errors in zero rank cases, and can potentionally give slightly
            ## faster and/or more accurate results:
            XXtYYt <- XXtYYt %*% XXtYYt

            ## Initial values:
            t.a.old <- Y0[,1]
            nit <- 0                        # for debugging
            repeat {
                nit <- nit + 1              # for debugging
                t.a <- XXtYYt %*% t.a.old
                t.a <- t.a / sqrt(c(crossprod(t.a)))
                if (sum(abs((t.a - t.a.old) / t.a), na.rm = TRUE) < tol)
                    break
                else
                    t.a.old <- t.a
                if (nit >= maxit) {         # for debugging
                    warning("No convergence in ", maxit, " iterations\n")
                    break
                }
            }
            nits[a] <- nit                  # for debugging

            u.a <- YYt %*% t.a
            utmp <- u.a / c(crossprod(t.a, u.a))
            wpw <- sqrt(c(crossprod(utmp, XXt) %*% utmp))
            TT[,a] <- t.a * wpw
            U[,a] <- utmp * wpw

            G <- In - tcrossprod(t.a)
            XXt <- G %*% XXt %*% G
            YYt <- G %*% YYt %*% G

            if (!stripped) Xresvar[a] <- sum(diag(XXt))
        }

        W <- crossprod(X0, U)
        W <- W / rep(sqrt(colSums(W * W)), each = npred)

        TTtTinv <- TT %*% diag(1 / colSums(TT * TT), ncol = ncol(TT))
        P <- crossprod(X0, TTtTinv)
        Q <- crossprod(Y0, TTtTinv)

        ## EM Step 2: Update estimates
        for (ind in seq_len(nrow(Xmissing))) {
            row <- Xmissing[ind,1]
            col <- Xmissing[ind,2]
            X[row,col] <- sum(TT[row,] * P[col,]) + Xmeans[col]
        }
        ## FIXME: make it a function, then do Y as well

        ## EM: Check for convergence
        changeX[EM.nit] <- mean(abs((X[Xmissing] - Xold) / Xold), na.rm = TRUE)
        change[EM.nit] <-sum(abs((abs(P) - abs(P.old)) / P), na.rm = TRUE)
        Xs[[EM.nit]] <- Xold
        Ps[[EM.nit]] <- P
        if (changeX[EM.nit] < EM.tol) {
            change <- change[1:EM.nit]
            changeX <- changeX[1:EM.nit]
            Xs[[EM.nit + 1]] <- X[Xmissing]
            break
        } else
            P.old <- P
        if (EM.nit >= EM.maxit) {       # for debugging
            warning("No EM convergence in ", EM.maxit, " iterations\n")
            break
        }
    }

    ## Calculate rotation matrix:
    if (ncomp == 1) {
        ## For 1 component, R == W:
        R <- W
    } else {
        PW <- crossprod(P, W)
        ## It is known that P^tW is right bi-diagonal (one response) or upper
        ## triangular (multiple responses), with all diagonal elements equal to 1.
        if (nresp == 1) {
            ## For single-response models, direct calculation of (P^tW)^-1 is
            ## simple, and faster than using backsolve.
            PWinv <- diag(ncomp)
            bidiag <- - PW[row(PW) == col(PW)-1]
            for (a in 1:(ncomp - 1))
                PWinv[a,(a+1):ncomp] <- cumprod(bidiag[a:(ncomp-1)])
        } else {
            PWinv <- backsolve(PW, diag(ncomp))
        }
        R <- W %*% PWinv
    }

    ## Calculate regression coefficients:
    for (a in 1:ncomp) {
        B[,,a] <- tcrossprod(R[,1:a, drop=FALSE], Q[,1:a, drop=FALSE])
    }

    if (stripped) {
        ## Return as quickly as possible
        list(coefficients = B, Xmeans = Xmeans, Ymeans = Ymeans)
    } else {
        ## Fitted values, residuals etc:
        for (a in 1:ncomp)
            fitted[,,a] <- tcrossprod(TT[,1:a, drop=FALSE], Q[,1:a, drop=FALSE])
        residuals <- - fitted + c(Y0)
        fitted <- fitted + rep(Ymeans, each = nobj) # Add mean
        Xvar <- diff(-c(Xtotvar, Xresvar))

        ## Add dimnames:
        objnames <- dnX[[1]]
        if (is.null(objnames)) objnames <- dnY[[1]]
        prednames <- dnX[[2]]
        respnames <- dnY[[2]]
        compnames <- paste("Comp", 1:ncomp)
        nCompnames <- paste(1:ncomp, "comps")
        dimnames(TT) <- dimnames(U) <- list(objnames, compnames)
        dimnames(R) <- dimnames(W) <- dimnames(P) <-
            list(prednames, compnames)
        dimnames(Q) <- list(respnames, compnames)
        dimnames(B) <- list(prednames, respnames, nCompnames)
        dimnames(fitted) <- dimnames(residuals) <-
            list(objnames, respnames, nCompnames)
        names(Xvar) <- compnames
        class(TT) <- class(U) <- "scores"
        class(P) <- class(W) <- class(Q) <- "loadings"

        list(coefficients = B,
             scores = TT, loadings = P,
             loading.weights = W,
             Yscores = U, Yloadings = Q,
             projection = R,
             Xmeans = Xmeans, Ymeans = Ymeans,
             fitted.values = fitted, residuals = residuals,
             Xvar = Xvar, Xtotvar = Xtotvar,
             nits = nits, EM.nit = EM.nit, change = change, Ps = Ps,
             changeX = changeX, Xs = Xs) # for debugging
    }
}
