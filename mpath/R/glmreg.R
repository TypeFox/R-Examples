#.First.lib <- function(lib, pkg)
.onLoad <- function(lib, pkg)
{
  library.dynam("mpath", pkg, lib)
}
                                        #coordinate descent algorithm
                                        #ref: Regularization paths for generalized linear models via coordinate descent, Friedman et al., JSS, 2010, 33(1)
glmreg <- function(x, ...) UseMethod("glmreg")

glmreg.default <- function(x, ...) {
    if (extends(class(x), "Matrix"))
        return(glmreg.matrix(x = x, ...))
    stop("no method for objects of class ", sQuote(class(x)),
         " implemented")
}

glmreg.formula <- function(formula, data, weights, offset=NULL, contrasts=NULL, 
  x.keep=FALSE, y.keep=TRUE, ...){
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

  RET <- glmreg_fit(X[,-1], Y, weights,...)
  RET$call <- match.call()
  RET <- c(RET, list(formula=formula, terms = mt, data=data,
   contrasts = attr(X, "contrasts"),
                       xlevels = .getXlevels(mt, mf)))
  if(x.keep) RET$x <- data[,colnames(data)%in%attr(mt, "term.labels")]
  if(y.keep) RET$y <- Y
  class(RET) <- "glmreg"
  RET
}
glmreg.matrix <- function(x, y, weights, offset=NULL, ...){
  RET <- glmreg_fit(x, y, weights,...)
  RET$call <- match.call()
  return(RET)
}

glmreg_fit <- function(x, y, weights, start=NULL, etastart=NULL, mustart=NULL, nlambda=100, lambda=NULL, lambda.min.ratio=ifelse(nobs<nvars,.05, .001),alpha=1, gamma=3, rescale=TRUE, standardize=TRUE, penalty.factor = rep(1, nvars),thresh=1e-6, eps.bino=1e-5, maxit=1000, eps=.Machine$double.eps, theta, family=c("gaussian", "binomial", "poisson", "negbin"), penalty=c("enet","mnet","snet"), convex=FALSE, x.keep=FALSE, y.keep=TRUE, trace=FALSE){
  if(!is.null(start) && !is.null(etastart) && !is.null(mustart))
  stop("start, etastart and mustart is for testing only\n")
  family <- match.arg(family)
  penalty <- match.arg(penalty)
  if(!family %in% c("gaussian", "binomial", "poisson", "negbin")){
    print(family)
    stop("'family' not recognizied\n")
  }
  if(family == "gaussian") rescale <- FALSE
  #if(family!="gaussian" && y < 0)
  #  stop("except for gaussian family, response should be nonnegative")
    if (gamma <= 1 && penalty=="mnet") stop("gamma must be greater than 1 for the mnet penalty")
    if (gamma <= 2 && penalty=="snet") stop("gamma must be greater than 2 for the snet penalty")
    if (alpha < 0 || alpha > 1){
    cat("alpha=", alpha)
  stop("alpha must be greater than 0 and less than or equal to 1")
}
    if (alpha == 0 && is.null(lambda)){
  stop("not designed for alpha=0 and lambda=NULL\n")
}
  if(!is.null(etastart) && !is.null(mustart)){
   if((length(etastart) != length(mustart)) || length(etastart) != length(y))
  stop("length of etastart and mustart should be the same as y\n") 
}
  ## avoid any problems with 1D or nx1 arrays by as.vector.
  nm <- dim(x)
  nobs <- n <- nm[1]
  nvars <- m <- nm[2]
if(missing(weights)) weights=rep(1,nobs)
  weights <- as.vector(weights)
  penalty <- match.arg(penalty)
  if(!is.null(weights) && !is.numeric(weights))
    stop("'weights' must be a numeric vector")
  ## check weights and offset
  if( !is.null(weights) && any(weights < 0) ){
    stop("negative weights not allowed")
}
  if(family=="binomial"){
   if(is.factor(y))
   y <- as.integer(y) - 1
   if(any(y > 1))
    stop("response should be between 0 and 1 for family=binomial")
  }
  if(family=="negbin"){
  if(missing(theta))
  stop("theta has to be provided for family='negbin'\n")
  else if(length(theta) > 1)
  stop("theta has to be one scale parameter for family='negbin'\n")
  else if(theta <= 0){
  cat("theta=", theta, "\n")
  stop("theta has to be positive for family='negbin'\n")
}
}
  if(family != "negbin") theta <- 1 ### this theta is not useful but required as input in Fortran glmlink subroutine
  tracetype <- 0
  if(trace)
    tracetype <- 1
 famtype <- switch(family,
      "gaussian"=1,
      "binomial"=2,
      "poisson"=3,
      "negbin"=4)
 pentype <- switch(penalty,
      "enet"=1,
      "mnet"=2,
      "snet"=3)
  stantype <- as.numeric(standardize)
                                        #  warning("alp = 1/theta in glm negative.binomial function\n")
  if(length(penalty.factor) != nvars) stop("number of penalty.factor should be the same as nvars")
  im <- inactive <- seq(m)
### compute the pathwise coordinate descent, cf: section 2.5 in Friedman et al., JSS 2010, 33(1)
  if(is.null(weights)) weights <- rep(1, n)
  wt <- weights/sum(weights)
  if(is.null(mustart) || is.null(etastart)){
  tmp <- init(wt, y, family=family)
  mu <- tmp$mu
  eta <- rep(tmp$eta,n)
}
  else{
  mu <- mustart
  eta <- etastart
  }
  if(is.null(lambda)){
  tmp <- init(wt, y, family=family)
  mu <- tmp$mu
  eta <- rep(tmp$eta,n)
    w <- .Fortran("glmlink",
      n=as.integer(1),
      mu=as.double(mu),
      family=as.integer(famtype),
      theta=as.double(theta),
      w = as.double(0),
      ep = as.double(eps.bino),
      package="mpath")$w
    z <- .Fortran("zeval",
      n=as.integer(n),
      y=as.double(y),
      eta=as.double(eta),
      mu=as.double(rep(mu,n)),
      w=as.double(rep(w,n)),
      family=as.integer(famtype),
      z=as.double(rep(0,n)),
      package="mpath")$z
    lmax <- findlam(x=x, y=y, weights=weights, family=family, theta=theta, mu=mu, w=w, z=z, alpha=alpha, penalty.factor=penalty.factor, standardize=standardize) 
#    if(penalty %in% c("mnet", "snet") && !rescale) lmax <- 0.5 * sqrt(lmax)
    lpath <- seq(log(lmax), log(lambda.min.ratio * lmax), length.out=nlambda)
    lambda <- exp(lpath)
  }
  else nlambda <- length(lambda)
### ref: McCullagh and Nelder, 2nd edition, 1989, page 121
  if(length(mu) != n)
  mu <- rep(mu, n)
  nulldev <- .Fortran("deveval",
  n=as.integer(n),
  y=as.double(y),
  mu=as.double(mu),
  theta=as.double(theta),
  weights=as.double(weights),
  family=as.integer(famtype),
  dev=as.double(0),
  package="mpath")$dev
  beta <- matrix(0, ncol=nlambda, nrow=m)
  if(is.null(start))
  startv <- 0
  if(!is.null(start)){
   if(length(start) != (m+1))
   stop("length of start doesn't match x dimension\n")
   else{
   startv <- 1
   if(length(start) > 1){
    for(j in 1: nlambda)
    beta[,j] <- start[-1]
}
  b0 <- rep(start[1], nlambda)
}
}
  else{
  b0 <- rep(0, nlambda)
  start <- rep(0, dim(x)[2]+1)
}
  resdev <- rep(0, nlambda)
  yhat <- matrix(0, nobs, nlambda)
  penfac <- penalty.factor/sum(penalty.factor) * nvars
  lam <- penfac %o% lambda
  if(family=="gaussian") {
   innermaxit <- maxit
         maxit <- 1
   }
   else innermaxit <- 1
 wtnew <- weights/sum(weights)
  meanx <- drop(wtnew %*% x)
  if(standardize){
    xx <- scale(x, meanx, FALSE) # centers x
    xx <- sqrt(wtnew) * xx
    one <- rep(1,n)  
    normx <- sqrt(one %*% xx^2)
    xx <- scale(x, meanx, normx)
  }
  tmp <- .Fortran("outloop",
   x=as.double(x),
   y=as.double(y),
   weights=as.double(weights),
   wt=as.double(wt),
   n=as.integer(n),
   m=as.integer(m),
   penalty=as.integer(pentype),
   nlambda=as.integer(nlambda),
   lam=as.double(lam),
   alpha=as.double(alpha),
   gam=as.double(gamma),
   theta=as.double(theta),
   rescale=as.integer(rescale),
   mu=as.double(mu),
   eta=as.double(eta),
   family=as.integer(famtype),
   standardize=as.integer(stantype),
   nulldev=as.double(nulldev),
   thresh=as.double(thresh),
   maxit=as.integer(maxit),
   innermaxit=as.integer(innermaxit),
   eps=as.double(eps),
   trace=as.integer(tracetype),
   start=as.double(start),
   startv=as.integer(startv),
   b=as.double(beta),
   bz=as.double(b0),
   resdev=as.double(resdev),
   ypre=as.double(yhat),
   convout=as.integer(rep(0,nlambda)), 
   satu=as.integer(0),
   good=as.integer(nlambda),
   ep = as.double(eps.bino),
   outpll = as.double(matrix(0, maxit, nlambda)),
   package="mpath") 
   if(nlambda > 1 || tmp$satu==0)
   good <- 1:tmp$good # only those non-saturated models
   else if(tmp$satu==1)
   return(RET <- list(satu=1))
### Names
  if(is.null(colnames(x))) varnames <- paste("V", 1:ncol(x), sep="")
  else varnames <- colnames(x)
  beta <- matrix(tmp$b, ncol=nlambda)[,good]
  b0 <- tmp$bz[good]
  ### note: pll was from standardized beta values if standardize=TRUE
  if(trace)
  pll <- matrix(tmp$outpll, ncol=nlambda)
  else pll <- NULL
  convex.min <- if (convex && all(weights==1)) convexMin(rbind(b0, beta), xx, penalty, gamma, lambda*(1-alpha), family, theta=theta) else NULL
  if(convex && any(weights!=1)) warnings("code for convex region implemented for weights=1 only\n")
   if (standardize){ 
   #if (family != "gaussian" && standardize){ 
  beta <- as.matrix(beta/as.vector(normx))
  b0 <- b0 - crossprod(meanx,beta)
  if (family == "gaussian")
   b0 <- mean(y) + b0    ### changed 4/22/2015
  }
  else normx <- NULL
  resdev <- tmp$resdev[good]
  yhat <- matrix(tmp$ypre, ncol=nlambda)[,good] 
  theta <- rep(theta, nlambda)
  if(is.null(dim(beta)))
  names(beta) <- colnames(x)
  else{
  rownames(beta) <- colnames(x)
  # colnames(beta) <- lambda
  colnames(beta) <- round(lambda,digits=4)
}
  RET <- list(family=family,standardize=standardize, satu=tmp$satu, lambda=lambda[good], nlambda=length(lambda[good]), beta=beta, b0=b0, meanx=meanx, normx=normx, theta=theta[good], nulldev=nulldev, resdev=resdev, pll = pll, fitted.values=yhat, converged=tmp$convout[good], convex.min=convex.min, penalty.factor=penalty.factor, gamma=gamma, alpha=alpha)
  if(x.keep) RET$x <- x
  if(y.keep) RET$y <- y
  class(RET) <- "glmreg"
  RET$twologlik <- try(2*logLik(RET, newx=x, y=y, weights=weights))
###penalized log-likelihood function value for rescaled beta
  penval <- .Fortran("penGLM", 
  start=as.double(beta),
  n=as.integer(n),
  m=as.integer(m),
  lambda=as.double(rep(RET$lambda, m)),
  alpha=as.double(alpha),
  gam=as.double(gamma),
  penalty=as.integer(pentype),
  pen=as.double(0.0),
  package="mpath")$pen
  RET$penval <- penval
  if(standardize)
  RET$pllres <- RET$twologlik/2 - n*penval
  else RET$pllres <- RET$twologlik/2 - penval
  if(nlambda == 1){
  RET$df <- sum(abs(beta) > 0) ### df= number of nonzero coefficients (intercept excluded)
  RET$aic <- -RET$twologlik + 2*(1+RET$df)
  RET$bic <- -RET$twologlik + log(n)*(1+RET$df)
}
  else{
  RET$df <- apply(abs(beta) > 0, 2, sum) ##number of nonzero coefficients for each lambda
  RET$aic <- -RET$twologlik + 2*(1+RET$df) #intercept included
  RET$bic <- -RET$twologlik + log(n)*(1+RET$df) #intercept included
}
  RET  
}


g <- function(mu, family, eps.bino=1e-5){
  switch(family,
         "gaussian"={
           eta <- mu
         },
         "binomial"={
           ep <- eps.bino
           n <- length(mu)
           eta <- rep(0, n)
           for(i in 1:n){
           if(1-mu[i] > ep && mu[i] > ep)
           eta[i] <- log(mu[i]/(1-mu[i]))
           }           
         },
         "poisson"={
           eta <- log(mu)
         },
         "negbin"={
           eta <- log(mu)
         }
        )
 return(eta)
}

### eta is the estimated beta_0 in the intercept-only model
init <- function(wt, y, family){
  mu <- wt %*% y
  switch(family,
         "gaussian"={
           eta <- mu
         },
         "binomial"={
           eta <- log(mu/(1-mu))
         },
         "poisson"={
           eta <- log(pmax(1, mu))
         },
         "negbin"={ ###Negative binomial regression by Hilbe, page 52, Table 4.1. Note, Table 8.3 has errors.
           eta <- log(pmax(1, mu))
         }
         )
  list(mu=mu, eta=eta)
}

findlam <- function(x, y, weights, family, theta=NULL,mu, w, z, alpha, penalty.factor, eps.bino=1e-5, standardize=TRUE){ 
  if (alpha <= 0 || alpha > 1){
  stop("alpha must be greater than 0 and less than or equal to 1, but alpha=", alpha)
}
  if(standardize)
  weights <- weights/sum(weights)
  penfac <- penalty.factor/sum(penalty.factor) * dim(x)[2]
  acvar <- which(penfac > 0)
### regression from unregularied variables
  if(length(acvar) < length(penalty.factor)){
    x <- x[, acvar]
    penfac <- penfac[acvar]
  }

  if(family == "negbin" && is.null(theta)){
  fit <- glm.nb(y ~ 1, weights = weights)
  mu <- fit$fitted.values
  w <- .Fortran("glmlink",
  n = as.integer(length(y)),
  mu = as.double(mu),
  family = as.integer(4),
  theta = as.double(fit$theta),
  w = as.double(rep(0, length(y))),
  ep = as.double(eps.bino),
  package="mpath")$w
  resid <- .Fortran("zeval",
  n = as.integer(length(y)),
  y = as.double(y),
  eta = as.double(rep(0, length(y))),
  mu = as.double(mu),
  w = as.double(w),
  family = as.integer(4),
  z = as.double(rep(0, length(y))),
  package="mpath")$z
}
  else{
  if(length(acvar) < length(penalty.factor))
    resid <- lm(z ~ x[,-acvar], weights=weights)$resid
  else resid <- z-weighted.mean(z, weights)
}
  if(standardize){
    xx <- x
     meanx <- drop(weights %*% x)
    x <- scale(x, meanx, FALSE) # centers x
    x <- sqrt(weights) * x
    one <- rep(1, length(y))
    normx <- sqrt(drop(one %*% (x^2)))
    x <- scale(xx, meanx, normx)
    lmax <- max(abs(crossprod(x,weights*w*resid))/(penfac*alpha))
     }
  else
    lmax <- max(abs(crossprod(x,weights*w*resid))/(penfac*alpha))
  lmax
}

deviance.glmreg <- function(object, newx, y, weights, na.action=na.pass, ...){
if(missing(newx) && missing(y)) return(object$resdev)
family <- object$family
if(missing(weights)) weights <- rep(1, length(y))
if(!is.null(object$terms)){
 mf <- model.frame(delete.response(object$terms), newx, na.action = na.action, xlev = object$xlevels)
 newx <- model.matrix(delete.response(object$terms), mf, contrasts = object$contrasts)
 newx <- newx[,-1] ### remove the intercept
}
object$terms <- NULL ### now newx is matrix, but predict function requires a dataframe for newx unless terms is NULL
mu <- predict(object, newx, type="response")
 famtype <- switch(family,
      "gaussian"=1,
      "binomial"=2,
      "poisson"=3,
      "negbin"=4)
dev <- rep(NA, object$nlambda)
for(j in 1:object$nlambda)
dev[j] <- .Fortran("deveval",
n=as.integer(length(y)),
y=as.double(y),
mu=as.double(mu[,j]),
theta=as.double(object$theta[j]),
weights=as.double(weights),
family= as.integer(famtype),
dev=as.double(0),
package="mpath")$dev
return(dev)
}

### convert glm object to class glmreg, thus can be used to predict newdata
                conv2glmreg <- function(object, family=c("poisson", "negbin")){
                family <- match.arg(family)
                if(class(object)=="glm") class(object) <- "glmreg"
                object$family <- family
                object$nlambda <- 1
                namecount <- names(object$coefficients)
                object$coefficients <- matrix(object$coefficients, ncol=1)
                rownames(object$coefficients) <- namecount
                object
}

