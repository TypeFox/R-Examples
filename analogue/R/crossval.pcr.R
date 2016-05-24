`crossval.pcr` <- function(obj, method = c("LOO","kfold","bootstrap"),
                           ncomp, nboot = 100, nfold = 10, folds = 5,
                           verbose = getOption("verbose"), ...) {
    method <- match.arg(method)
    x <- obj$data$x
    y <- obj$data$y
    N <- NROW(x)
    M <- NCOL(x)

    FUN <- obj$tranFun
    PARMS <- obj$tranParms

    ##B <- coef(obj)

    if(identical(method, "LOO")) {
        nr <- N-1 ## number of rows - 1 for LOO
        ## form ncomp, as LOO we have potentially 1 less component than usual
        ncomp <- if(missing(ncomp)) {
            min(nr - 1, M) ## uses nr which already has 1 removed
        } else {
            if(ncomp < 1 || ncomp > (newcomp <- min(nr - 1, M))) {
                warning("'ncomp' inappropriate for LOO CV. Resetting to max possible.")
                newcomp
            } else {
                ncomp
            }
        }
        ## matrix of predictions
        pred <- matrix(ncol = ncomp, nrow = N)
        ## display progress
        if(verbose) {
            writeLines("\n  LOO Cross-validation:")
            pb <- txtProgressBar(min = 0, max = nr, style = 3)
            on.exit(close(pb))
            on.exit(cat("\n"), add = TRUE)
        }
        nc <- seq_len(ncomp)
        for(i in seq_len(N)) {
            if(verbose)
                setTxtProgressBar(pb, i)
            ## apply transformation to X[-i, ]
            TRAN <- obj$tranFun(x[-i, , drop = FALSE])
            X <- TRAN$data
            ## apply transformation to X[i, ], using parms from above
            Xi <- obj$tranFun(x[i, , drop = FALSE], apply = TRUE,
                              parms = TRAN$parms)$data
            ## centre the training data
            Xbar <- colMeans(X)
            ybar <- mean(y[-i])
            X <- sweep(X, 2, Xbar, "-")
            Y <- y[-i] - ybar
            ## fit model to subset
            FIT <- fitPCR(X = X, Y = Y, ncomp = ncomp, n = nr, m = M)
            ## predict for 1:ncomps components
            for(j in nc) {
                B0 <- obj$yMean - obj$xMeans %*% FIT$B[, j, drop = FALSE]
                pred[i, j] <- Xi %*% FIT$B[, j, drop = FALSE] + B0
            }
        }
    }
    if(identical(method, "kfold")) {
        ## form ncomp, as k-fold we have ceiling(N / nfold) fewer sites
        maxComp <- min((N - ceiling(N / nfold)) - 1, M)
        ncomp <- if(missing(ncomp)) {
            maxComp## uses nr which already has 1 removed
        } else {
            if(ncomp < 1 || ncomp > maxComp) {
                warning("'ncomp' inappropriate for k-fold CV. Resetting to max possible.")
                maxComp
            } else {
                ncomp
            }
        }
        pred <- array(NA, dim = c(N, ncomp, folds))
        if(verbose) {
            writeLines("\n   n k-fold Cross-validation:")
            ii <- 1
            pb <- txtProgressBar(min = 0, max = folds * nfold, style = 3)
            on.exit(close(pb))
            on.exit(cat("\n"), add = TRUE)
        }
        ind <- rep(seq_len(nfold), length = N) ## k-fold group indicator
        nc <- seq_len(ncomp)
        ## this is the n in n k-fold CV, allowing n repeated k-folds
        for(i in seq_len(folds)) {
            ## do a k-fold CV
            pind <- ind[sample.int(N, N, replace = FALSE)]
            ## the main k-fold CV loop
            for(k in seq_len(nfold)) {
                if(verbose) {
                    setTxtProgressBar(pb, ii)
                    ii <- ii + 1
                }
                sel <- pind == k   ## sel is samples in leave out group
                N.oob <- sum(sel) ## N in leave out group
                N.mod <- sum(!sel)  ## N in the model
                sel <- which(sel) ## convert to indices
                ## apply transformation to X[-sel, ]
                TRAN <- obj$tranFun(x[-sel, , drop = FALSE])
                X <- TRAN$data
                ## apply transformation to X[sel, ], using parms from above
                Xi <- obj$tranFun(x[sel, , drop = FALSE], apply = TRUE,
                                  parms = TRAN$parms)$data
                ## centre the training data
                Xbar <- colMeans(X)
                ybar <- mean(y[-sel])
                X <- sweep(X, 2, Xbar, "-")
                Y <- y[-sel] - ybar
                ## fit model to subset
                FIT <- fitPCR(X = X, Y = Y, ncomp = ncomp, n = N.mod, m = M)
                ## predict for 1:ncomps components
                for(j in nc) {
                    B0 <- obj$yMean - obj$xMeans %*% FIT$B[, j, drop = FALSE]
                    pred[sel, j, i] <- Xi %*% FIT$B[, j, drop = FALSE] + rep(B0, N.oob)
                }
            }
        }
    }
    if(identical(method, "bootstrap")) {
        ## form ncomp, as if this was a standard training set setting.
        maxComp <- min(N - 1, M)
        ncomp <- if(missing(ncomp)) {
            maxComp
        } else {
            if(ncomp < 1 || ncomp > maxComp) {
                warning("'ncomp' inappropriate for bootstrap CV. Resetting to max possible.")
                maxComp
            } else {
                ncomp
            }
        }
        pred <- array(NA, dim = c(N, ncomp, nboot))
        if(verbose) {
            writeLines("\n   Bootstrap Cross-validation:")
            pb <- txtProgressBar(min = 0, max = nboot, style = 3)
            on.exit(close(pb))
            on.exit(cat("\n"), add = TRUE)
        }
        ind <- seq_len(N) ## indicator for samples
        nc <- seq_len(ncomp)
        for(i in seq_len(nboot)) {
            if(verbose)
                setTxtProgressBar(pb, i)
            bSamp <- sample.int(N, N, replace = TRUE)
            sel <- which(!ind %in% bSamp) ## need indices!!!
            N.oob <- NROW(x[sel, , drop = FALSE])
            ##N.mod <- N - N.oob ## not sure I need this
            ## apply transformation to X[-sel, ]
            TRAN <- obj$tranFun(x[bSamp, , drop = FALSE])
            X <- TRAN$data
            ## apply transformation to OOB samples, using parms from above
            Xi <- obj$tranFun(x[sel, , drop = FALSE], apply = TRUE,
                              parms = TRAN$parms)$data
            ## centre the training data
            Xbar <- colMeans(X)
            ybar <- mean(y[bSamp])
            X <- sweep(X, 2, Xbar, "-")
            Y <- y[bSamp] - ybar
            ## fit model to subset
            FIT <- fitPCR(X = X, Y = Y, ncomp = ncomp, n = N, m = M)
            ## predict for 1:ncomps components
            for(j in nc) {
                B0 <- obj$yMean - obj$xMeans %*% FIT$B[, j, drop = FALSE]
                pred[sel, j, i] <- Xi %*% FIT$B[, j, drop = FALSE] +
                    rep(B0, N.oob)
            }
        }
    }

    ## fitted values and derived stats
    if(identical(method, "none")) {
        fitted <- pred
    } else if(method %in% c("kfold","bootstrap")) {
        fitted <- rowMeans(pred, na.rm = TRUE, dims = 2)
    } else {
        fitted <- rowMeans(pred)
    }
    residuals <- y - fitted                               ## residuals
    maxBias <- apply(residuals, 2, maxBias, y, n = 10)    ## maxBias
    avgBias <- colMeans(residuals)                        ## avgBias
    r2 <- apply(fitted, 2, cor, y)
    ## s1 & s2 components for model and training set
    ns <- rowSums(!is.na(pred), dims = 2)
    s1.train <- sqrt(rowSums((pred - as.vector(fitted))^2,
                             na.rm = TRUE, dims = 2) / as.vector(ns - 1))
    s1 <- sqrt(colMeans(s1.train^2))
    s2 <- sqrt(colMeans(residuals^2, na.rm = TRUE))
    ## s2.train <- sweep(pred, 1, y, "-")
    ## s2.train <- sqrt(rowMeans(s2.train^2, na.rm = TRUE, dims = 2))

    ## RMSEP
    ## rmsep.train <- sqrt(s1.train^2 + s2.train^2)
    rmsep2 <- sqrt(s1^2 + s2^2)
    rmsep <- sqrt(colMeans(residuals^2, na.rm = TRUE))
    fill <- rep(NA, ncomp)
    performance <- data.frame(comp = seq_len(ncomp),
                              R2 = r2,
                              avgBias = avgBias,
                              maxBias = maxBias,
                              RMSEP = rmsep,
                              RMSEP2 = fill,
                              s1 = fill,
                              s2 = fill)
    ## add in the bits we can only do if bootstrapping or multiple fold
    ## k-fold CV
    if(identical(method, "bootstrap") ||
       (identical(method, "kfold") && folds > 1)) {
        performance$s1 <- s1
        performance$s2 <- s2
        performance$RMSEP2 <- rmsep2
    }

    ## more additions to the call
    .call <- match.call()
    .call[[1]] <- as.name("crossval")

    ## return object
    out <- list(fitted.values = fitted,
                residuals = residuals,
                ##rmsep = rmsep.train, ## technically not in crossval.wa
                ##s1 = s1.train,       ## so shoould they be in here?
                ##s2 = s2.train,       ## need to formalise the crossval class
                performance = performance,
                call = .call,
                CVparams = list(method = method, nboot = nboot,
                nfold = nfold, folds = folds))

    class(out) <- "crossval"
    out
}
