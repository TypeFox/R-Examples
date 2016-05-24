## crossval method for wa()
`crossval.wa` <- function(obj, method = c("LOO","kfold","bootstrap"),
                          nboot = 100, nfold = 10, folds = 5,
                          verbose = getOption("verbose"), ...) {
    method <- match.arg(method)
    X <- obj$orig.x
    ENV <- obj$orig.env
    N <- NROW(X)
    M <- NCOL(X)
    tolOpts <- obj$options.tol
    Dtype <- obj$deshrink
    if(identical(method, "LOO")) {
        pred <- numeric(N)
        nr <- N-1 ## number of rows - 1 for LOO
        if(verbose) {
            writeLines("\n  LOO Cross-validation:")
            pb <- txtProgressBar(min = 0, max = nr, style = 3)
            on.exit(close(pb))
            on.exit(cat("\n"), add = TRUE)
        }
        for(i in seq_along(pred)) {
            if(verbose)
                setTxtProgressBar(pb, i)
            opt <- w.avg(X[-i, ], ENV[-i])
            if(obj$tol.dw)
                pred[i] <- predWAT(X, ENV, i, opt, tolOpts, nr, M,
                                   Dtype)
            else
                pred[i] <- predWA(X, ENV, i, opt, Dtype)
        }
    }
    if(identical(method, "kfold")) {
        oob.pred <- matrix(NA, ncol = folds, nrow = N)
        if(verbose) {
            writeLines("\n   n k-fold Cross-validation:")
            pb <- txtProgressBar(min = 0, max = folds, style = 3)
            on.exit(close(pb))
            on.exit(cat("\n"), add = TRUE)
        }
        ind <- rep(seq_len(nfold), length = N) ## k-fold group indicator
        ## n k-fold s
        for(i in seq_len(folds)) {
            if(verbose)
                setTxtProgressBar(pb, i)
            ## do a k-fold CV
            pind <- ind[sample.int(N, N, replace = FALSE)] ## sure this should be replace = FALSE
            for(k in seq_len(nfold)) {
                sel <- pind == k   ## sel is samples in leave out group
                N.oob <- sum(sel) ## N in leave out group
                N.mod <- sum(!sel)  ## N in the model
                sel <- which(sel) ## convert to indices
                opt <- w.avg(X[-sel, , drop = FALSE], ENV[-sel])
                if(obj$tol.dw) {
                    oob.pred[sel, i] <- predWAT(X, ENV, sel, opt, tolOpts,
                                                N.mod, M, Dtype)
                } else {
                    oob.pred[sel, i] <- predWA(X, ENV, sel, opt, Dtype)
                }
            }
        }
        pred <- rowMeans(oob.pred, na.rm = TRUE)
    }
    if(identical(method, "bootstrap")) {
        oob.pred <- matrix(NA, ncol = nboot, nrow = N)
        if(verbose) {
            writeLines("\n   Bootstrap Cross-validation:")
            pb <- txtProgressBar(min = 0, max = nboot, style = 3)
            on.exit(close(pb))
            on.exit(cat("\n"), add = TRUE)
        }
        ind <- seq_len(N) ## indicator for samples
        for(i in seq_len(nboot)) {
            if(verbose)
                setTxtProgressBar(pb, i)
            bSamp <- sample.int(N, N, replace = TRUE)
            sel <- which(!ind %in% bSamp) ## need indices!!!
            N.oob <- NROW(X[sel, , drop = FALSE])
            N.mod <- N - N.oob
            opt <- w.avg(X[-sel, , drop = FALSE], ENV[-sel])
            if(obj$tol.dw)
                oob.pred[sel, i] <- predWAT(X, ENV, sel, opt, tolOpts,
                                             N.mod, M, Dtype)
            else
                oob.pred[sel, i] <- predWA(X, ENV, sel, opt, Dtype)
        }
        pred <- rowMeans(oob.pred, na.rm = TRUE)
    }
    resid <- ENV - pred
    out <- list(fitted.values = pred, residuals = resid)
    performance <- data.frame(R2 = cor(pred, ENV)^2,
                              avgBias = mean(resid),
                              maxBias = unname(maxBias(resid, ENV)),
                              RMSEP = sqrt(mean(resid^2)),
                              RMSEP2 = NA,
                              s1 = NA,
                              s2 = NA)
    if(identical(method, "bootstrap") ||
       (identical(method, "kfold") && folds > 1)) {
        performance$s1 <- sqrt(mean(apply(oob.pred, 1, sd, na.rm = TRUE)^2))
        performance$s2 <- sqrt(mean(resid^2))
        performance$RMSEP2 <- sqrt(performance$s1^2 + performance$s2^2)
    }
    out$performance <- performance
    .call <- match.call()
    .call[[1]] <- as.name("crossval")
    out$call <- .call
    out$CVparams <- list(method = method, nboot = nboot, nfold = nfold,
                         folds = folds)
    class(out) <- "crossval"
    out
}

`predWAT` <- function(X, ENV, i, optima, tolOpts, nr, nc,
                      deSh) {
    tol <- w.tol(X[-i, ], ENV[-i], optima, tolOpts$useN2)
    tol <- fixUpTol(tol, tolOpts$na.tol, tolOpts$small.tol,
                    tolOpts$min.tol, tolOpts$f, ENV[-i])
    wa.env <- WATpred(X[-i, ], optima, tol, nr, nc)
    p <- WATpred(X[i, , drop = FALSE], optima, tol, 1, nc)
    deMod <- deshrink(ENV[-i], wa.env, deSh)
    deshrinkPred(p, coef(deMod), deSh)
}

`predWA` <- function(X, ENV, i, optima, deSh) {
    wa.env <- WApred(X[-i, ], optima)
    p <- WApred(X[i, , drop = FALSE], optima)
    deMod <- deshrink(ENV[-i], wa.env, deSh)
    deshrinkPred(p, coef(deMod), deSh)
}
