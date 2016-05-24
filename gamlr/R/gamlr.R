##########################################################
##### path estimation for log penalized regression  ######
##########################################################

## Wrapper function; most happens in c
gamlr <- function(x, y, 
            family=c("gaussian","binomial","poisson"),
            gamma=0, 
            nlambda=100, 
            lambda.start=Inf,  
            lambda.min.ratio=0.01, 
            free=NULL, 
            standardize=TRUE, 
            obsweight=NULL,
            varweight=NULL,
            prexx=(p<500),  
            tol=1e-7, 
            maxit=1e5,
            verb=FALSE, ...)
{
  on.exit(.C("gamlr_cleanup", PACKAGE = "gamlr"))

  ## integer family codes
  family=match.arg(family)
  famid = switch(family, 
    "gaussian"=1, "binomial"=2, "poisson"=3)

  ## data checking (more follows)
  y <- checky(y,family)
  n <- length(y)

  # observation weights
  if(!is.null(obsweight))
    if(family!="gaussian"){ 
      warning("non-null obsweight are ignored for family!=gaussian")
      obsweight <- NULL }
  if(is.null(obsweight)) obsweight <- rep(1,n)
  stopifnot(all(obsweight>0))
  stopifnot(length(obsweight)==n)

  ## extras
  xtr = list(...)

  ## aliases from glmnet or other gamlr terminology
  if(!is.null(xtr$thresh)) tol = xtr$thresh
  if(!is.null(xtr$lmr)) lambda.min.ratio = xtr$lmr
  if(!is.null(xtr$scale)) standardize = xtr$scale
  if(!is.null(xtr[['fix']])) xtr$shift = xtr$fix

  ## max re-weights
  if(is.null(xtr$maxrw)) xtr$maxrw = maxit # practically inf
  maxrw = xtr$maxrw

  ## fixed shifts 
  eta <- rep(0.0,n)
  if(!is.null(xtr$shift)){
    if(family=="gaussian") y = y-xtr$shift
    else eta <- xtr$shift   } 
  stopifnot(length(eta)==n)
  eta <- as.double(eta)

  ## get x dimension and names
  if(is.null(x)){
    if(any(c(family!="gaussian",
      is.null(xtr$vxx),
      is.null(xtr$vxy),
      is.null(xtr$vxsum),
      is.null(xtr$xbar))))
        stop("xx,xy,xsum,xbar are NULL or family!=`gaussian'; 
          this is not allowed if x=NULL")
    p <- length(xtr$xbar)
    varnames <- names(xtr$xbar) 
    x <- Matrix(0)
  } else{
    if(inherits(x,"numeric")) x <- matrix(x)
    if(inherits(x,"data.frame")) x <- as.matrix(x)
    if(inherits(x,"simple_triplet_matrix"))
      x <- sparseMatrix(i=x$i,j=x$j,x=x$v,
              dims=dim(x),dimnames=dimnames(x))
    p <- ncol(x)
    varnames <- colnames(x)
    stopifnot(nrow(x)==n) 
  }
  if(is.null(varnames)) varnames <- paste(1:p)

  # fixedcost (undocumented: additional fixed l1 penalty)
  if(is.null(xtr$fixedcost))
    xtr$fixedcost <- 0
  fixedcost = xtr$fixedcost
  if(length(fixedcost)!=p){ 
      fixedcost <- rep(fixedcost[1],p) }

  ## unpenalized columns
  if(length(free)==0) free <- NULL
  if(!is.null(free)){
    if(inherits(free,"character")){
      free <- na.omit(match(free,varnames))
      if(length(free)==0) free <- NULL
      print(free)}
    if(any(free < 1) | any(free>p)) stop("bad free argument.") 
    if(length(free)==p){
      nlambda <- 1
      lambda.start <- 0
    }
  }
  
  ## variable (penalty) weights
  if(is.null(varweight)) varweight <- rep(1,p)
  stopifnot(all(varweight>=0))
  stopifnot(length(varweight)==p)
  varweight[free] <- 0

  ## check and clean all arguments
  stopifnot(lambda.min.ratio<=1)
  stopifnot(all(c(nlambda,lambda.min.ratio)>0))
  stopifnot(all(c(lambda.start)>=0))
  stopifnot(all(c(tol,maxit)>0))
  if(lambda.start==0){
    nlambda <- 1
    standardize <- 0 }
  lambda <- double(nlambda)
  lambda[1] <- lambda.start

  ## stepsize
  delta <- exp( log(lambda.min.ratio)/(nlambda-1) )

  ## adaptation
  stopifnot(all(gamma>=0))
  if(length(gamma)==1){ 
    gamvec <- rep(gamma,p)
  } else{ 
    gamvec <- gamma }
  stopifnot(length(gamvec)==p)
  gamvec[free] <- 0

  ## PREXX stuff
  prexx = (prexx | !is.null(xtr$vxx)) & (family=="gaussian")
  if(prexx){
    if(is.null(xtr$xbar))
      xtr$xbar <- colMeans(x)
    xbar <- as.double(xtr$xbar)
    if(is.null(xtr$vxsum))
      xtr$vxsum <- colSums(x*obsweight)
    vxsum <- as.double(xtr$vxsum)
    if(is.null(xtr$vxx))
      xtr$vxx <- crossprod(x*sqrt(obsweight))
    vxx <- Matrix(xtr$vxx,sparse=FALSE,doDiag=FALSE)
    vxx <- as(vxx,"dspMatrix")
    stopifnot(ncol(vxx)==p)
    if(vxx@uplo=="L") vxx <- t(vxx)
    vxx <- as.double(vxx@x)
    if(is.null(xtr$vxy))
      xtr$vxy <- drop(crossprod(x,y*obsweight))
    vxy <- as.double(drop(xtr$vxy))
    stopifnot(all(p==
      c(length(xbar),length(vxsum),length(vxy))))
  } else{
    xbar <- double(p)
    vxsum <- double(p)
    vxy <- double(p)
    vxx <- double(0)
  }

  ## final x formatting
  x=as(x,"dgCMatrix") 
  stopifnot(all(is.finite(x@x)))


  ## drop it like it's hot
  fit <- .C("gamlr",
            famid=as.integer(famid), 
            n=n,
            p=p,
            l=length(x@i),
            xi=x@i,
            xp=x@p,
            xv=as.double(x@x),
            y=y,
            prexx=as.integer(prexx),
            xbar=xbar,
            xsum=vxsum,
            xx=vxx,
            xy=vxy,
            eta=eta,
            varweight=as.double(varweight),
            obsweight=as.double(obsweight),
            standardize=as.integer(standardize>0),
            nlambda=as.integer(nlambda),
            delta=as.double(delta),
            gamma=gamvec,
            fixedcost=as.double(fixedcost),
            tol=as.double(tol),
            maxit=as.integer(rep(maxit,nlambda)),
            maxrw=as.integer(rep(maxrw,nlambda)),
            lambda=as.double(lambda),
            deviance=double(nlambda),
            df=double(nlambda),
            alpha=as.double(rep(0,nlambda)),
            beta=as.double(rep(0,nlambda*p)),
            exits=integer(nlambda), 
            verb=as.integer(verb>0),
            PACKAGE="gamlr",
            NAOK=TRUE,
            dup=FALSE)

  ## coefficients
  nlambda <- fit$nlambda
  if(nlambda == 0){
    warning("could not converge for any lambda.")
    nlambda <- 1
    fit$deviance <- NA
  }

  if(any(fit$exits==1)) 
    warning("numerically perfect fit for some observations.")

  alpha <- head(fit$alpha,nlambda)
  names(alpha) <- paste0('seg',(1:nlambda))
  beta <- Matrix(head(fit$beta,nlambda*p),
                    nrow=p, ncol=nlambda, 
                    dimnames=list(varnames,names(alpha)),
                    sparse=TRUE)

  ## path stats
  lambda <- head(fit$lambda,nlambda)
  dev <- head(fit$deviance,nlambda)
  df <- head(fit$df,nlambda)
  names(df) <- names(dev) <- names(lambda) <- names(alpha)

  ## iterations
  iter <- cbind(ncycle=head(fit$maxit,nlambda),
                nreweight=head(fit$maxrw,nlambda))
  rownames(iter) <- names(alpha)

  ## build return object and exit
  out <- list(lambda=lambda, 
             gamma=gamma,
             nobs=fit$n,
             family=family,
             alpha=alpha,
             beta=beta, 
             df=df,
             deviance=dev,
             iter=iter,
             free=free,
             call=sys.call(1)) 

  class(out) <- "gamlr"
  invisible(out)
}

## just response argument checking
checky <- function(y, family){ 
  if(is.data.frame(y)) y <- as.matrix(y)
  y <- drop(y)
  stopifnot(is.null(dim(y)))
  if(is.factor(y)&family=="binomial") y <- as.numeric(y)-1
  stopifnot(all(is.finite(y)))
  return(as.double(y))
}


#### S3 method functions

plot.gamlr <- function(x, against=c("pen","dev"), 
                      col=NULL, 
                      select=TRUE, df=TRUE, ...)
{
  nlambda <- ncol(x$beta)
  if(!is.null(x$free)) 
    x$beta <- x$beta[-x$free,,drop=FALSE]
  p <- nrow(x$beta) 
  nzr <- unique(x$beta@i)+1
  if(length(nzr)==0 | (x$lambda[1]==0)) return("nothing to plot")
  beta <- as.matrix(x$beta[nzr,,drop=FALSE])

  if(!is.null(col)){
    if(length(col)==1) col <- rep(col,p)
    col <- col[nzr] }
  else col <- 1:6 # matplot default

  against=match.arg(against)
  if(against=="pen"){
      xv <- log(x$lambda)
      xvn <- "log lambda"
  } else if(against=="dev"){
      xv <- x$dev
      xvn <- "deviance"
  } else
    stop("unrecognized 'against' argument.")

  if(!is.finite(xv[1])) 
    stop("refusing to plot an unconverged fit")

  argl = list(...)
  if(is.null(argl$ylim)) argl$ylim=range(beta)
  if(is.null(argl$ylab)) argl$ylab="coefficient"
  if(is.null(argl$xlab)) argl$xlab=xvn
  if(is.null(argl$lty)) argl$lty=1
  if(is.null(argl$bty)) argl$bty="n"
  do.call(plot, c(list(x=xv, y=rep(0,nlambda), col="grey70", type="l"), argl))

  matplot(xv, t(beta), col=col, add=TRUE, type="l", lty=argl$lty)

  if(df){
    dfi <- unique(round(
      seq(1,nlambda,length=ceiling(length(axTicks(1))))))
    axis(3,at=xv[dfi], labels=round(x$df[dfi],1),tick=FALSE, line=-.5) }

  if(select){
    abline(v=xv[which.min(AICc(x))], lty=2, col="grey20")
  }
}

coef.gamlr <- function(object, select=NULL, k=2, ...)
{
  if(length(select)==0){
    select <- which.min(AICc(object,k=k))
    if(length(select)==0) select <- 1
  }
  else if(select==0)
   select <- 1:ncol(object$beta)

  select[select>ncol(object$beta)] <- ncol(object$beta)
  return(rBind(intercept=object$alpha[select], 
              object$beta[,select,drop=FALSE]))
}

predict.gamlr <- function(object, newdata,
                    type = c("link", "response"), ...)
{
  if(inherits(newdata,"data.frame")) 
    newdata <- as.matrix(newdata)
  if(inherits(newdata,"simple_triplet_matrix"))
    newdata <- sparseMatrix(
                    i=newdata$i,
                    j=newdata$j,
                    x=newdata$v,
                    dims=dim(newdata),
                    dimnames=dimnames(newdata))

  B <- coef(object, ...)
  eta <- matrix(B[1,],
            nrow=nrow(newdata),
            ncol=ncol(B),
            byrow=TRUE)
  eta <- eta + newdata%*%B[-1,,drop=FALSE]
                              
  type=match.arg(type)
  if(object$family=="binomial" & type=="response")
    return(apply(eta,2,function(c) 1/(1+exp(-c))))
  if(object$family=="poisson" & type=="response")
    return(apply(eta,2,exp))

  return(eta)
}

summary.gamlr <- function(object, ...){
  print(object)

  return(data.frame(
    lambda=object$lambda,
    par=diff(object$b@p)+1,
    df=object$df,
    r2=1-object$dev/object$dev[1],
    aicc=AICc(object)))
}

print.gamlr <- function(x, ...){
  cat("\n")
  cat(sprintf(
    "%s gamlr with %d inputs and %d segments.", 
    x$family, nrow(x$beta), ncol(x$beta)))
  cat("\n\n")
}

logLik.gamlr <- function(object, ...){
  if(object$family=="gaussian"){
    ll <- -0.5*object$nobs*log(object$dev/object$nobs)
  } else{ ll <- -0.5*object$dev }
  attr(ll,"nobs") = object$nobs
  attr(ll,"df") = object$df
  class(ll) <- "logLik"
  ll
}




