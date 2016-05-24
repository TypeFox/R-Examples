cv.glmregNB <- function(formula, data, weights, lambda=NULL, 
nfolds=10, foldid, plot.it=TRUE, se=TRUE, n.cores=2, 
...){
    call <- match.call()
 mf <- Call <- match.call()
                                     m <- match(c("formula", "data", "subset", "weights", "na.action",
                                                "etastart", "mustart", "offset"), names(mf), 0)
                                            mf <- mf[c(1, m)]
                                            mf$drop.unused.levels <- TRUE
                                            mf[[1L]] <- as.name("model.frame")
                                            mf <- eval.parent(mf)
                                            Terms <- attr(mf, "terms")
                                            Y <- model.response(mf, "numeric")
  ## null model support
                                            X <- if (!is.empty.model(Terms)) model.matrix(Terms, mf, contrasts) else matrix(,NROW(Y),0)
  nobs <- n <- length(Y)
  nvars <- m <- dim(X) - 1
 weights <- model.weights(mf)
 if(!length(weights)) weights <- rep(1, nrow(mf))
  if(any(weights < 0)) stop("negative weights not allowed")

  K <- nfolds
      glmregNB.obj <- do.call("glmregNB", list(formula, data, weights, lambda=lambda, ...))
    lambda <- glmregNB.obj$lambda
    nlambda <- length(lambda)
    if(missing(foldid))
    all.folds <- cv.folds(n, K)
    else all.folds <- foldid 
    fraction <- seq(nlambda)
    registerDoParallel(cores=n.cores)
    i <- 1  ###needed to pass R CMD check with parallel code below
    residmat <- foreach(i=seq(K), .combine=cbind) %dopar% {
      omit <- all.folds[[i]]
### changed 5/20/2013 fixed theta
      fitcv <- do.call("glmregNB", list(formula, data[-omit,], weights[-omit], nlambda=nlambda, lambda=lambda, theta.est=FALSE, theta0=glmregNB.obj$theta, ...))
### remove the first column, which is for intercept
      fitcv$terms <- NULL ### logLik requires data frame if terms is not NULL
      logLik(fitcv, newx=X[omit,-1, drop=FALSE], Y[omit], weights=weights[omit])
   }
    cv <- apply(residmat, 1, mean)
    cv.error <- sqrt(apply(residmat, 1, var)/K)
    lambda.which <- which.max(cv)
    obj<-list(fit=glmregNB.obj, residmat=residmat, fraction = fraction, cv = cv, cv.error = cv.error, foldid=all.folds, lambda.which= lambda.which, lambda.optim = lambda[lambda.which])
    class(obj) <- "cv.glmreg"
    if(plot.it) plot(obj,se=se)
    obj
  }

