`predict.pcr` <- function(object, newdata, ncomp = object$ncomp,
                          CV = c("none", "LOO", "bootstrap", "kfold"),
                          verbose = FALSE, nboot = 100, kfold = 10,
                          folds = 5, ...) {
    takeData <- function(x, new) {
        want <- (spp.names <- colnames(x$data$x)) %in% colnames(new)
        want <- spp.names[want]
        new[, want, drop = FALSE]
    }
    if(missing(newdata))
        return(fitted(object))
    ## store names of new samples
    newSamp <- rownames(newdata)
    newdata <- as.matrix(newdata)
    if (missing(CV))
        CV <- "none"
    CV <- match.arg(CV)
    Nnew <- NROW(newdata)
    N <- nrow(object$data$x)
    M <- ncol(object$data$x)

    ind <- seq_len(N) ## indicator for samples

    ## extract the training data & transformation fun
    trainX <- object$data$x
    trainY <- object$data$y
    tranFun <- object$tranFun

    if(identical(CV, "none")) {
        B <- coef(object)
        newdata <- takeData(object, newdata)
        ## apply transformation to newdata
        tf <- object$tranFun(newdata)
        newdata <- tf$data
        ## do predictions
        ## matrix of predictions
        pred <- matrix(ncol = ncomp, nrow = Nnew)
        for(j in seq_len(ncomp)) {
            B0 <- object$yMean - object$xMeans %*% B[, j]
            pred[, j] <- newdata %*% B[, j] + rep(B0, Nnew)
        }
    } else if (identical(CV, "LOO")) {
        nr <- N-1 ## for LOO
        M <- ncol(object$data$x)
        ## form ncomp, as LOO we have potentially 1 less component than usual
        ncomp <- if(missing(ncomp)) {
            min(nr - 1, M) ## uses nr which already has 1 removed
        } else {
            if(ncomp < 1 || ncomp > (newcomp <- min(nr - 1, M))) {
                warning("'ncomp' inappropriate for LOO CV.
Resetting to max possible.")
                newcomp
            } else {
                ncomp
            }
        }

        pred <- array(NA, dim = c(Nnew, ncomp, N))
        comps <- seq_len(ncomp)

        for (i in seq_len(N)) {
            ## which samples are in the training set
            loo <- ind[-i]

            ## do the heavy lifting
            pred[, , i] <-
                pcrCVfit(trainX, trainY, tf = tranFun, newdata,
                         parray = pred[, , i], take = loo,
                         N = nr, Nnew = Nnew, M = M, ncomp = ncomp,
                         comps = comps)
        }
        fitted <- rowMeans(pred, na.rm = TRUE, dims = 2)
        ## other computations on `pred` needed for SEs etc
        pred <- fitted
    } else if (identical(CV, "bootstrap")) {
        ## form ncomp, as if this was a standard training set setting.
        maxComp <- min(N - 1, M)
        ncomp <- if(missing(ncomp)) {
            maxComp
        } else {
            if(ncomp < 1 || ncomp > maxComp) {
                warning("'ncomp' inappropriate for bootstrap CV.
Resetting to max possible.")
                maxComp
            } else {
                ncomp
            }
        }

        pred <- array(NA, dim = c(Nnew, ncomp, nboot))
        comps <- seq_len(ncomp)

        for (i in seq_len(nboot)) {
            bSamp <- sample.int(N, N, replace = TRUE)

            ## do the heavy lifting
            pred[, , i] <-
                pcrCVfit(trainX, trainY, tf = tranFun, newdata,
                         parray = pred[, , i], take = bSamp,
                         N = N, Nnew = Nnew, M = M, ncomp = ncomp,
                         comps = comps)
        }
        fitted <- rowMeans(pred, na.rm = TRUE, dims = 2)
        ## other computations on `pred` needed for SEs etc
        pred <- fitted
    } else if (identical(CV, "kfold")) {
        ## form ncomp, as k-fold we have ceiling(N / nfold) fewer sites
        maxComp <- min((N - ceiling(N / kfold)) - 1, M)
        ncomp <- if(missing(ncomp)) {
            maxComp## uses nr which already has 1 removed
        } else {
            if(ncomp < 1 || ncomp > maxComp) {
                warning("'ncomp' inappropriate for k-fold CV.
Resetting to max possible.")
                maxComp
            } else {
                ncomp
            }
        }

        comps <- seq_len(ncomp)

        ## array for prediction
        pred <- array(NA, dim = c(Nnew, ncomp, folds, kfold))
        ind <- rep(seq_len(kfold), length = N) ## k-fold group indicator

        ## this is the n in n k-fold CV, allowing n repeated k-folds
        ##ii <- 0
        for(i in seq_len(folds)) {
            ## do a k-fold CV
            pind <- ind[sample.int(N, N, replace = FALSE)]
            ## the main k-fold CV loop
            for(k in seq_len(kfold)) {
                ##ii <- ii + 1
                kSamp <- pind != k
                Nk <- sum(kSamp) ## number of samples in training set
                kSamp <- which(kSamp)

                ## do the heavy lifting
                pred[, , i, k] <-
                    pcrCVfit(trainX, trainY, tf = tranFun,
                             newdata, parray = pred[, , i, k],
                             take = kSamp, N = Nk, Nnew = Nnew,
                             M = M, ncomp = ncomp, comps = comps)
            }
        }
        fitted <- rowMeans(pred, na.rm = TRUE, dims = 2)
        ## other computations on `pred` needed for SEs etc
        pred <- fitted
    }else {
        stop("Unknown crossvalidation method.")
    }
    rownames(pred) <- newSamp
    colnames(pred) <- paste0("PC", seq_len(ncomp))
    pred
}

## take is a vector of indices for samples to include in the training
## set during prediction
`pcrCVfit` <- function(X, Y, tf, newdata, parray, take, N, Nnew, M,
                       ncomp, comps) {
    ## apply transformation to training data
    TRAN <- tf(X[take, , drop = FALSE])
    X <- TRAN$data
    ## apply transformation to newdata, using parms from above
    Xnew <- tf(newdata, apply = TRUE, parms = TRAN$parms)$data
    ## centre the training data
    Xbar <- colMeans(X)
    ybar <- mean(Y[take])
    X <- sweep(X, 2, Xbar, "-")
    Y <- Y[take] - ybar
    ## fit model to subset
    FIT <- fitPCR(X = X, Y = Y, ncomp = ncomp, n = N, m = M)
    for(j in comps) {
        B0 <- ybar - Xbar %*% FIT$B[, j]
        parray[ , j] <- Xnew %*% FIT$B[, j] + rep(B0, Nnew)
    }
    parray
}
