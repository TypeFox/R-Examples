## Principal Components Regression Models

## generic
`pcr` <- function(x, ...) {
    UseMethod("pcr")
}

## default
`pcr.default` <- function(x, y, ncomp, tranFun, ...) {
    ## convert to matrices for speed
    x <- data.matrix(x)

    ## Store data
    XX <- x
    YY <- y

    ## dimensions
    Nx <- NROW(x)
    Ny <- length(y)
    Mx <- NCOL(x) ## number of predictor vars

    ## Apply a transformation if specified
    fun.supplied <- FALSE
    if(!missing(tranFun)) {
        tranFun <- match.fun(tranFun)
        x <- tranFun(x)
        tranParms <- x$parms
        x <- x$data
        attr(x, "parms") <- NULL
        fun.supplied <- TRUE
    } else {
        tranFun <- function(x, ...) list(data = x, parms = NA)
        tranParms <- list(NA)
        ##tranFun <- "none"
    }

    ## centre x and y
    xMeans <- colMeans(x)
    yMean <- mean(y)
    x <- sweep(x, 2, xMeans, "-")
    y <- y - yMean

    ## How many components?
    if(missing(ncomp)) {
        ncomp <- min(Nx - 1, Mx)
    } else {
        if(ncomp < 1 || ncomp > (newcomp <- min(Nx - 1, Mx))) {
            warning("Invalid 'ncomp'. Resetting to max possible.")
            ncomp <- newcomp
        }
    }

    ## fit via wrapper
    FIT <- fitPCR(X = x, Y = y, ncomp = ncomp, n = Nx, m = Mx)

    ## other model output
    residuals <- y - FIT$Yhat
    fitted.values <- sweep(FIT$Yhat, 2, yMean, "+")

    ## model performance
    Y <- y + yMean
    performance <- data.frame(R2 = drop(cor(fitted.values, Y)),
                              avgBias = colMeans(residuals),
                              maxBias = apply(residuals, 2, maxBias, Y),
                              RMSE = sqrt(colMeans(residuals^2)))

    ## apply some names to prettify to output
    rownames(performance) <- colnames(FIT$B) <- colnames(fitted.values) <-
        colnames(residuals) <- names(FIT$varExpl) <- colnames(FIT$TT) <-
            colnames(FIT$P) <- rownames(FIT$tQ) <- paste0("PC", seq_len(ncomp))
    rownames(FIT$B) <- rownames(FIT$P) <- colnames(x)
    rownames(fitted.values) <- rownames(residuals) <- rownames(FIT$TT) <- rownames(x)

    ## get and fix up the call
    .call <- match.call()
    .call[[1]] <- as.name("pcr")

    ## return object
    Obj <- list(fitted.values = fitted.values,
                coefficients = FIT$B,
                residuals = residuals,
                scores = FIT$TT,
                loadings = FIT$P,
                Yloadings = t(FIT$tQ),
                xMeans = xMeans,
                yMean = yMean,
                varExpl = FIT$varExpl,
                totvar = sum(x * x),
                call = .call,
                tranFun = tranFun,
                tranParms = tranParms,
                performance = performance,
                ncomp = ncomp,
                data = list(x = XX, y = YY))
    class(Obj) <- "pcr"
    Obj
}

`pcr.formula` <- function(formula, data, subset, na.action,
                          ..., model = FALSE) {
    ## the function call
    .call <- match.call()
    ## need to reset due to method dispatch
    .call[[1]] <- as.name("pcr")

    ## keep only the arguments which should go into the model frame
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1]] <- as.name("model.frame")
    mf <- eval.parent(mf)
    attr(attr(mf, "terms"), "intercept") <- 0

    ## terms objects
    mt <- attr(mf, "terms")
    ## model matrices
    y <- model.response(mf, "numeric")
    x <- model.matrix(mt, mf)

    ## fit & the PCR
    Obj <- pcr.default(x = x, y = y, ...)
    Obj$na.action <- attr(mf, "na.action")
    Obj$terms <- mt
    if(model)
        Obj$model <- mf
    Obj$call <- .call
    Obj
}

`Hellinger` <- function(x, ...) {
    list(data = tran(x, method = "hellinger"), parms = list(NA))
    ## ignore 'apply' as no meta-parameters
}

`ChiSquare` <- function(x, apply = FALSE, parms) {
    if(apply) {
        ## apply pre-computed meta-parameters to transform
        ## test samples to match training samples
        ## take only variables in x for which we have parms
        ## match on attr(parms, "variables")
        data <- sqrt(parms$gsum) * x / outer(parms$rsum, sqrt(parms$csum))
        res <- list(data = data, parms = parms)
    } else {
        ## perform transformation and preserve the meta-parameters
        if (any(x < 0, na.rm = TRUE)) {
            k <- min(x, na.rm = TRUE)
            warning("'x'contains negative entries: result may be nonsense\n")
        } else {
            k <- .Machine$double.eps
        }
        gsum <- sum(x, na.rm = TRUE)
        rsum <- pmax(k, rowSums(x, na.rm = TRUE))
        csum <- colSums(x, na.rm = TRUE)
        parms <- list(gsum = gsum, rsum = rsum, csum = csum)
        attr(parms, "variables") <- colnames(x)
        res <- list(data = sqrt(gsum) * x/outer(rsum, sqrt(csum)),
                    parms = parms)
    }
    res ##return
}

`print.pcr` <- function(x, digits = min(getOption("digits"), 4), ...) {
    cat("\n")
    writeLines(strwrap("Principal Component Regression Model", prefix = "\t"),
               sep = "\n\n")
    writeLines(strwrap("Call:"))
    print(x$call)
    cat("\n")
    writeLines(strwrap(paste("No. of Components:", x$ncomp)), sep = "\n\n")
    writeLines("RMSE (Apparent):")
    perf <- x$performance[, "RMSE"]
    names(perf) <- rownames(x$performance)
    print(perf)
}
