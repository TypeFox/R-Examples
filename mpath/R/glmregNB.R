###Adapted from file MASS/R/negbin.R
glmregNB <- function(formula, data, weights, nlambda=100, lambda=NULL, lambda.min.ratio=ifelse(nobs<nvars,.05, .001),alpha=1, gamma=3, rescale=TRUE, standardize=TRUE, penalty.factor = rep(1, nvars), thresh=1e-3, maxit.theta=25, maxit=1000, eps=.Machine$double.eps, trace=FALSE, start = NULL, etastart = NULL, mustart = NULL, theta.est=TRUE, theta0=NULL, init.theta=ifelse(theta.est, theta0[1],NULL), link=log, penalty=c("enet","mnet","snet"), method="glmreg_fit", model=TRUE, x.keep=FALSE, y.keep=TRUE, contrasts=NULL, convex=FALSE, ...)
{
  if(!theta.est)
   if(length(theta0)!=nlambda)
   stop("length of theta0 must be the same as nlambda if theta.est=FALSE\n")
   else if(any(theta0 <= 0))
   stop("theta0 must be positive\n") 
   penalty <- match.arg(penalty)
  if(penalty=="enet") convex=FALSE
  loglik <- function(n, th, mu, y, w)
    sum(w*(lgamma(th + y) - lgamma(th) - lgamma(y + 1) + th * log(th) +
           y * log(mu + (y == 0)) - (th + y) * log(th + mu)))
    link <- substitute(link)
    fam0 <- if(missing(init.theta))
        do.call("poisson", list(link = link))
    else
        do.call("negative.binomial", list(theta = init.theta, link = link))
   dots <- list(...)
  mf <- Call <- match.call()
                                        m <- match(c("formula", "data", "subset", "weights", "na.action",
                                                "etastart", "mustart", "offset"), names(mf), 0)
                                            mf <- mf[c(1, m)]
                                            mf$drop.unused.levels <- TRUE
                                            mf[[1L]] <- as.name("model.frame")
                                            mf <- eval.parent(mf)
                                            Terms <- attr(mf, "terms")
                                            if(method == "model.frame") return(mf)
                                            Y <- model.response(mf, "numeric")
  ## null model support
                                            X <- if (!is.empty.model(Terms)) model.matrix(Terms, mf, contrasts) else matrix(,NROW(Y),0)
                                        w <- model.weights(mf)
  nm <- dim(X[,-1])
  nobs <- n <- nm[1]
  nvars <- m <- nm[2]
  if(is.null(penalty.factor))
  penalty.factor <- rep(1, nvars)
#  if(missing(weights)) w <- weights <- rep(1, length(y))
  if(!length(w)) w <- weights <- rep(1, nrow(mf))
  else w <- weights
  if(any(w < 0)) stop("negative weights not allowed")
   offset <- model.offset(mf)
    ## these allow starting values to be expressed in terms of other vars.
    mustart <- model.extract(mf, "mustart")
    etastart <- model.extract(mf, "etastart")
  ### thoughts: use glm.nb to estimate intercept-only model which can take care of theta estimate 
  if(is.null(lambda)){
     lmax <- findlam(x=X[,-1], y=Y, weights=weights, family="negbin", theta=NULL, mu=NULL, w=NULL, z=NULL, alpha=alpha, penalty.factor=penalty.factor, standardize=standardize)
    lpath <- seq(log(lmax), log(lambda.min.ratio * lmax), length.out=nlambda)
    lambda <- exp(lpath)
  }
  else 
  nlambda <- length(lambda)
if(standardize){
    wt <- weights/sum(weights)
    one <- rep(1, length(Y))
    Xnew <- X[,-1]
    if(dim(Xnew)[2] > 0){
    xx <- Xnew
    meanx <- drop(wt %*% Xnew)
    Xnew <- scale(Xnew, meanx, FALSE) # centers x
    Xnew <- sqrt(wt) * Xnew
    normx <- sqrt(drop(one %*% (Xnew^2)))
    X[,-1] <- Xnew <- scale(xx, meanx, normx)
}
}
  beta <- matrix(NA, ncol=nlambda, nrow=m)
  fitted <- matrix(NA, ncol=nlambda, nrow=n)
  pll <- matrix(NA, ncol=nlambda, nrow=maxit)
  b0 <- tht <- nulldev <- resdev <- rep(NA, nlambda)                           
            #    dots <- list(...)
  k <- 1
  convout <- twologlik <- rep(NA, nlambda)
  while(k <= nlambda){  
    if(trace) message("loop in lambda:", k, "\n")
    if(k==1){
    if(trace){ cat("Initial fit family is:")
    print(fam0)
}
    if(trace) message("Initial fit:")
   if(!missing(method)) {
        if(!exists(method, mode = "function"))
            stop("unimplemented method: ", sQuote(method))
        glm.fitter <- get(method)
    } else {
        method <- "glm.fit"
        glm.fitter <- stats::glm.fit
    }
    if(is.null(mustart) || missing(init.theta)){
    fit <- glm.nb(Y ~ 1, weights=weights)
    mu <- fit$fitted.values
    th <- fit$theta
}
    else{
    mu <- mustart
    th <- init.theta
}
    if(trace)
      message("Initial value for theta:", signif(th))
}
    if(!theta.est) th <- theta0[k]
    iter <- 0
    d1 <- sqrt(2 * max(1, fit$df.residual))
    d2 <- del <- 1
    Lm <- loglik(n, th, mu, Y, w)
    converged <- FALSE                                    
    while((iter <- iter + 1) <= maxit.theta && !converged){
      eta <- log(mu)
       fit <- glmreg_fit(x=X[,-1], y=Y, weights=w, lambda=lambda[k],alpha=alpha,gamma=gamma,rescale=rescale, standardize=standardize, penalty.factor = penalty.factor, thresh=thresh, maxit=maxit, eps=eps, family="negbin", theta=th, trace=trace, penalty=penalty)
       #fit <- glmreg_fit(x=X[,-1], y=Y, weights=w, start = start, etastart=eta, mustart = mu, lambda=lambda[k],alpha=alpha,gamma=gamma,rescale=rescale, standardize=standardize, penalty.factor = penalty.factor, thresh=thresh, maxit=maxit, eps=eps, family="negbin", theta=th, trace=trace, penalty=penalty)
      t0 <- th
      mu <- fit$fitted.values
      if(theta.est){
#      th <- theta.ml(Y, mu, sum(w), w, limit=maxit.theta,
      th <- theta.ml(Y, mu, sum(w), w,
                     trace = trace)
#      if(inherits(th, "try-error")) browser()
}
      else th <- theta0[k]
      start <- c(fit$b0, fit$beta)
      del <- t0 - th
      Lm0 <- Lm
      penval <- ifelse(standardize, n*fit$penval, fit$penval)
      Lm <- loglik(n, th, mu, Y, w) - penval
      fit$df.residual <- n - fit$df - 1
      d1 <- sqrt(2 * max(1, fit$df.residual))
      converged <- abs((Lm0 - Lm)/d1) + abs(del/d2) < 1e-8
      if(trace) {
        Ls <- loglik(n, th, Y, Y, w)
        Dev <- 2 * (Ls - Lm)
        message("Theta(", iter, ") =", signif(th),
                ", 2(Ls - Lm) =", signif(Dev))
      }
    }
    if(!is.null(attr(th, "warn"))) fit$th.warn <- attr(th, "warn")
    if(trace && iter > maxit.theta) {
      warning("alternation limit reached")
      fit$th.warn <- gettext("alternation limit reached")
    }
    tht[k] <- th
    beta[,k] <- as.vector(fit$beta)
    b0[k] <- fit$b0
    Lm <- loglik(n, th, mu, Y, w)
    twologlik[k] <- as.vector(2 * Lm)
    nulldev[k] <- fit$nulldev
    resdev[k] <- fit$resdev
    if(trace) pll[,k] <- fit$pll
    convout[k] <- converged
    fitted[,k] <- fit$fitted.values
    k <- k + 1
  }
  class(fit) <- c("glmregNB", "glmreg", "lm")
  fit$terms <- Terms
  fit$formula <- as.vector(attr(Terms, "formula"))
  ## make result somewhat reproducible
  Call$init.theta <- signif(as.vector(th), 10)
    Call$link <- link
    fit$call <- Call
    if(model) fit$model <- mf
    fit$na.action <- attr(mf, "na.action")
    if(x.keep) fit$x <- data[,colnames(data)%in%attr(Terms, "term.labels")]
    if(!y.keep) fit$y <- NULL
 if(standardize) {
  beta <- beta/normx
  b0 <- b0 - crossprod(meanx, beta)
}
    fit$twologlik <- twologlik
  fit$df <- apply(abs(beta) > 0, 2, sum) ### excluding intercept and theta
  fit$aic <- -fit$twologlik + 2*(2+fit$df)
  fit$bic <- -fit$twologlik + log(n)*(2+fit$df)
  fit$fitted.values <- fitted
    fit$contrasts <- attr(X, "contrasts")
    fit$xlevels <- .getXlevels(Terms, mf)
    fit$method <- method
    fit$offset <- offset
 if(is.null(colnames(X[,-1]))) varnames <- paste("V", 1:ncol(X[,-1]), sep="")
  else varnames <- colnames(X[,-1])
 dimnames(beta) <- list(varnames, round(lambda, digits=4))
    fit$beta <- beta
    fit$b0 <- matrix(b0, nrow=1)

    fit$lambda <- lambda
    fit$nlambda <- nlambda
    fit$nulldev <- nulldev
    fit$resdev <- resdev
    if(trace)
    fit$pll <- pll
    else fit$pll <- NULL
    fit$theta <- as.vector(tht)
    fit$converged <- convout
    fit
}


summary.glmregNB <- function(object, ...)
{
  if(!is.null(object$th.warn))
  summ <- object[c("theta", "twologlik", "th.warn")]
  else summ <- object[c("theta", "twologlik")]
  class(summ) <- "summary.glmregNB"
  summ
}

print.summary.glmregNB <- function(x, ...)
{
  NextMethod()
  cat("\n              Theta: ", format(round(x$theta, 3)),
      "\n")
  if(!is.null(x$th.warn))
    cat("Warning while fitting theta:", x$th.warn,"\n")
  cat("\n 2 x log-likelihood: ", format(round(x$twologlik, 3)), "\n")
  invisible(x)
}
pen_sum <- function(b, n, lambda, alpha, gamma, penalty){
lc <- length(b)
pen1 <- rep(0, lc)
if(lc > 0){
for(j in 1:lc)
pen1[j] <- pen_eval(abs(b[j]), lone=lambda*alpha, ltwo=lambda*(1-alpha), gamma=gamma, penalty=penalty)
}
return(n*(sum(pen1)))
}

