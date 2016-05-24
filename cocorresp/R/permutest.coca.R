"permutest.coca" <- function(x, R0 = NULL, permutations = 99,
                             n.axes = x$n.axes, verbose = TRUE, ...) {
    permtest <- function(Y, X1, X0 = NULL, permutations, step) {
        Y.dim <- dim(Y)
        X1.dim <- dim(X1)
        if(!identical(Y.dim[1], X1.dim[1])) {
            stop("Matrix X1 must have the same number of rows as Y")
        }
        if (is.null(X0)) {
            mu <- matrix(0, nrow = Y.dim[1], ncol = Y.dim[2])
        } else {
            X0.dim <- dim(X0)
            if(!identical(Y.dim[1], X0.dim[1]))
                stop("Matrix X0 must have the same number of rows as Y")
            if((X1.dim[2] <= X0.dim[2]))
                stop("Matrix X1 must have more columns than X0")
            mu <- qr.coef(qr(X0), Y)
            mu <- X0 %*% mu
        }
        E <- Y - mu
        SS0 <- sum(sum(E^2))
        stati <- numeric(length = permutations + 1)
        stati[1] <- teststat(Y = E, X0 = X0, X1 = X1, step)
        YresPerm <- matrix(0, nrow = Y.dim[1], ncol = Y.dim[2])
        for(i in 2:(permutations + 1)) {
            YresPerm[sample(Y.dim[1]), ] <- E
            stati[i] <- teststat(Y = YresPerm, X0 = X0, X1 = X1, step)
        }
        pval <- sum(stati >= stati[1]) / (permutations + 1)
        retval <- list(pval = pval, stat = stati[1], stati = stati)
        class(retval) <- "permtest"
        return(retval)
    }
    teststat <- function(Y, X0, X1, step) {
        Psi <- coinertiaI(X = Y, Y = X1, fast = TRUE)
        if(is.null(X0)) {
            fit.X0 <- 0
            X0X1 <- Psi[ , 1, drop = FALSE]
        } else {
            fit.X0 <- residualMatrix(Y = Y, X = X0)$inertia$fitted
            X0X1 <- matrix(c(X0, Psi[ , 1, drop = FALSE]), ncol = step)
        }
        resid.res <- residualMatrix(Y = Y, X = X0X1)
        fit.lambda1 <- (fit.X0X1 <- resid.res$inertia$fitted) - fit.X0
        retval <- fit.lambda1 / (resid.res$inertia$total - fit.X0X1)
        return(retval)
    }
    residualMatrix <- function(Y, X) {
        Q <- qr.coef(qr(X), Y)
        Yf <- X %*% Q
        Yr <- Y - Yf
        tot.inertia <- sum(sum(Y^2))
        resid.inertia <- sum(sum(Yr^2))
        fit.inertia <- sum(sum(Yf^2))
        zerosum <- (fit.inertia + resid.inertia) - tot.inertia
        if (abs(zerosum) > 0.000001)
            warning("Residual inertia + fitted inertia did not equal total inertia.\n\t",
                    call. = FALSE)
        retval <- list(inertia = list(total = tot.inertia,
                       residual = resid.inertia,
                       fitted = fit.inertia),
                       Yr = Yr, zerosum = zerosum)
        class(retval) <- "residualMatrix"
        return(retval)
    }
    if(!inherits(x, "predcoca"))
        stop("x must be of class 'predcoca'")
    if(is.null(R0)) {
        .R0 <- x$R0
    } else {
        .R0 <- R0
    }
    Ychi1 <- x$Ychi$Ychi1
    Ychi2 <- x$Ychi$Ychi2
    if(n.axes > x$n.axes) {
        n.axes <- x$n.axes
        warning("n.axes too large, reset to x$n.axes.")
    }
    pval <- permstat <- inertia <- fitax <- numeric(n.axes)
    for(j in 1:n.axes) {
        if(verbose) {
            cat("Permutations for axis:", j)
            flush.console()
        }
        if(j == 1)
            covar <- NULL
        ptest <- permtest(Ychi1, Ychi2, X0 = covar, permutations, step = j)
        if(j == 1)
            stati.Ax1 <- ptest$stati
        permstat[j] <- ptest$stati[1]
        pval[j] <- ptest$pval
        Psi <- coinertiaI(X = Ychi1, Y = Ychi2, fast = TRUE)[, 1, drop = FALSE]
        res.mat1 <- residualMatrix(Ychi1, Psi)
        Ychi1 <- res.mat1$Yr
        if(j == 1) {
            total.inertia1 <- res.mat1$inertia$total
        }
        res.mat2 <- residualMatrix(Ychi2, Psi)
        Ychi2 <- res.mat2$Yr
        if(is.null(covar)) {
            covar <- cbind(NULL, Psi)
        } else {
            covar <- cbind(covar, Psi)
        }
        inertia[j] <- res.mat1$inertia$total
        fitax[j] <- res.mat1$inertia$fitted
        if(verbose) {
            cat(" - completed\n")
            flush.console()
        }
    }
    pcent.fit <- 100 * fitax / total.inertia1
    retval <- list(pval = pval, permstat = permstat,
                   total.inertia = total.inertia1,
                   inertia = inertia, fitax = fitax,
                   pcent.fit = pcent.fit, n.axes = n.axes,
                   ##Ychi1 = Ychi1, Ychi2 = Ychi2, stati.Ax1 = stati.Ax1
                   call = match.call())
    class(retval) <- "permutest.coca"
    retval
}

