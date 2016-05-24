cv.zipath <- function(formula, data, weights, nlambda=100, lambda.count=NULL, lambda.zero=NULL, 
                      nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=2, 
                      ...){
    call <- match.call()
    if(missing(foldid) && nfolds < 3)
        stop("smallest nfolds should be 3\n")
    nm <- dim(data)
    nobs <- n <- nm[1]
    nvars <- m <- nm[2]
    mf <- Call <- match.call()
    m <- match(c("formula", "data", "subset", "weights", "na.action",
                 "offset"), names(mf), 0)
    mf <- mf[c(1, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    Y <- data[,all.vars(terms(formula))[[1]]]

    ## null model support
    weights <- model.weights(mf)
    if(!length(weights)) weights <- rep(1, length(Y))
    if(any(weights < 0)) stop("negative weights not allowed")

    K <- nfolds
    zipath.obj <- do.call("zipath", list(formula, data, weights, nlambda=nlambda, lambda.count=lambda.count, lambda.zero=lambda.zero, ...))
    lambda.count <- zipath.obj$lambda.count
    lambda.zero <- zipath.obj$lambda.zero
    nlambda <- zipath.obj$nlambda
    if(missing(foldid))
        all.folds <- cv.folds(n, K)
    else {
        all.folds <- foldid
        K <- nfolds <- length(foldid)
    }
    fraction <- seq(nlambda)
    bic <- matrix(NA, nlambda, K)
    registerDoParallel(cores=n.cores)
    i <- 1  ###needed to pass R CMD check with parallel code below
    residmat <- foreach(i=seq(K), .combine=cbind) %dopar% {
        omit <- all.folds[[i]]
        fitcv <- do.call("zipath", list(formula, data[-omit,], weights[-omit], lambda.count=lambda.count, lambda.zero=lambda.zero, nlambda=nlambda, ...))
        logLik(fitcv, newdata=data[omit,, drop=FALSE], Y[omit], weights=weights[omit])
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    lambda.which <- which.max(cv)
    obj<-list(fit=zipath.obj, residmat=residmat, bic=bic, fraction = fraction, cv = cv, cv.error = cv.error, foldid=all.folds, lambda.which= lambda.which, lambda.optim = list(count=lambda.count[lambda.which], zero=lambda.zero[lambda.which]))
    class(obj) <- c("cv.zipath", "cv.glmreg")
    if(plot.it) plot(obj,se=se)
    obj
}

coef.cv.zipath <- function(object, which=object$lambda.which, model = c("full", "count", "zero"), ...) {
    model <- match.arg(model)
    rval <- object$fit$coefficients
    rval <- switch(model,
                   "full" = list(count=rval$count[, which], zero=rval$zero[, which]),
                   "count" = rval$count[,which],
                   "zero" = rval$zero[,which])
    rval
}
