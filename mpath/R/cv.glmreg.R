"cv.folds" <-
    function(n, folds = 10)
    {
        split(sample(1:n), rep(1:folds, length = n))
    }

cv.glmreg <- function(x, ...) UseMethod("cv.glmreg", x)

cv.glmreg.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(cv.glmreg.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

cv.glmreg.formula <- function(formula, data, weights, offset=NULL, ...){
    ## extract x, y, etc from the model formula and frame
    if(!attr(terms(formula, data=data), "intercept"))
        stop("non-intercept model is not implemented")
    if(missing(data)) data <- environment(formula)
    mf <- match.call(expand.dots = FALSE)
    m <- match(c("formula", "data", "weights",
                 "offset"), names(mf), 0L)
    mf <- mf[c(1L, m)]
    mf$drop.unused.levels <- TRUE
    mf[[1L]] <- as.name("model.frame")
    mf <- eval(mf, parent.frame())

    mt <- attr(mf, "terms") # allow model.frame to have updated it

    Y <- model.response(mf, "any") # e.g. factors are allowed
    ## avoid problems with 1D arrays, but keep names
    if(length(dim(Y)) == 1L) {
        nm <- rownames(Y)
        dim(Y) <- NULL
        if(!is.null(nm)) names(Y) <- nm
        if(!is.null(nm)) names(Y) <- nm
    }
    ## null model support
    X <- if (!is.empty.model(mt)) model.matrix(mt, mf, contrasts) else matrix(,NROW(Y), 0L)
    ## avoid any problems with 1D or nx1 arrays by as.vector.
    weights <- as.vector(model.weights(mf))
    if(!length(weights)) weights <- rep(1, nrow(mf))
    else if(any(weights < 0)) stop("negative weights not allowed")
    if(!is.null(weights) && !is.numeric(weights))
        stop("'weights' must be a numeric vector")
    if(length(weights) != length(Y))
        stop("'weights' must be the same length as response variable")

    offset <- as.vector(model.offset(mf))
    if(!is.null(offset)) {
        if(length(offset) != NROW(Y))
            stop(gettextf("number of offsets is %d should equal %d (number of observations)", length(offset), NROW(Y)), domain = NA)
    }
### End of addition 08/07/2012 

    RET <- cv.glmreg_fit(X[,-1], Y, weights,...)
    RET$call <- match.call()
    return(RET)
}
cv.glmreg.matrix <- function(x, y, weights, offset=NULL, ...){
    RET <- cv.glmreg_fit(x, y, weights,...)
    RET$call <- match.call()
    return(RET)
}

cv.glmreg_fit <- function(x, y, weights, lambda=NULL, balance=TRUE, 
                          family=c("gaussian", "binomial", "poisson", "negbin"), 
                          nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=2, 
                          ...){
    call <- match.call()
    if(missing(foldid) && nfolds < 3)
        stop("smallest nfolds should be 3\n")
    family <- match.arg(family)
    nm <- dim(x)
    nobs <- n <- nm[1]
    nvars <- m <- nm[2]
    if(missing(weights)) weights <- rep(1, nobs)
    K <- nfolds
    glmreg.obj <- glmreg_fit(x, y, weights, lambda=lambda, family=family, ...)
    lambda <- glmreg.obj$lambda
    nlambda <- length(lambda)
    if(missing(foldid)){
        if(family=="binomial" && balance)  
            all.folds <- balanced.folds(y, K)
        else all.folds <- cv.folds(length(y), K)
    }
    else all.folds <- foldid
    fraction <- seq(nlambda)
    registerDoParallel(cores=n.cores)
    i <- 1  ###needed to pass R CMD check with parallel code below
    residmat <- foreach(i=seq(K), .combine=cbind) %dopar% {
        omit <- all.folds[[i]]
        fitcv <- glmreg_fit(x[ - omit,,drop=FALSE ], y[ -omit], weights=weights[- omit], lambda=lambda, family=family, ...)
	logLik(fitcv, newx=x[omit,, drop=FALSE], y[omit], weights=weights[omit])
    }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    lambda.which <- which.max(cv)
    obj<-list(fit=glmreg.obj, residmat=residmat, fraction = fraction, cv = cv, cv.error = cv.error, foldid=all.folds, lambda.which= lambda.which, lambda.optim = lambda[lambda.which])
    if(plot.it) plot.cv.glmreg(obj,se=se)
    class(obj) <- "cv.glmreg"
    obj
}

"plot.cv.glmreg" <-
    function(x,se=TRUE,ylab=NULL, main=NULL, width=0.02, col="darkgrey", ...){
        fraction <- x$fraction
        cv <- x$cv
        cv.error <- x$cv.error
        if(is.null(ylab))
            ylab <- "log-likelihood"
        plot(log(fraction), cv, type = "b", xlab = expression(log(lambda)), ylab= ylab, ylim = range(cv, cv + cv.error, cv - cv.error), main=main)
        if(se)
            error.bars(log(fraction), cv + cv.error, cv - cv.error,
                       width = width, col=col)

        invisible()
    }

coef.cv.glmreg=function(object,which=object$lambda.which,...){
    coef(object$fit,which=which,...)
}
