#############################################################
#	loggammarob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: February, 5, 2015
#	Version: 0.2
#	Copyright (C) 2015 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob <- function(x, start=NULL, weights=rep(1, length(x)), method=c("oneWL", "WQTau", "WL", "QTau", "ML"), control, ...) {
  method <- match.arg(method)
  x <- na.omit(x)
  if (!is.null(na <- attr(x, "na.action")))
    weights <- weights[-na]
  or <- order(x)
  if (any(or!=1:length(x))) {
    x <- sort(x)
    warning("data 'x' are sorted")
  }
  weights <- weights[or]
  
  if (missing(control)) 
    control <- if (missing(method))
                 loggammarob.control(...)
               else
                 loggammarob.control(method = method, ...)
    cl <- match.call()

  if (!missing(control) && !missing(method) && method != control$method) {
    warning("Methods argument set by method is different from method in control\n", "Using method = ", method)
    control$method <- method
  }
  ## x.orig <- x
  ## mx <- median(x)
  ## sx <- mad(x)
  ## x <- (x - mx)/sx
  
  if (control$method=="QTau")
    result <- loggammarob.QTau(x, weights, control)
  else if (method=="WQTau")
    result <- loggammarob.WQTau(x, weights, control)
  else if (method=="WL")
    result <- loggammarob.WL(x, start, weights, control)
  else if (method=="oneWL")
    result <- loggammarob.oneWL(x, start, weights, control)
  else if (method=="ML")
    result <- loggammarob.ML(x, start, weights, control)
  ## result$mu <- result$mu*sx+mx
  ## result$sigma <- result$sigma*sx  
  result$eta <- Exp.response(result$mu, result$sigma, result$lambda)
  ## result$data <- x.orig
  result$data <- x  
  result$method <- method
  result$call <- cl
  class(result) <- "loggammarob"
  return(result)
}

#############################################################
#	loggammarob.control function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.control <- function(method="oneWL", tuning.rho=1.547647, tuning.psi=6.08, lower=-7, upper=7, n=201, max.it=750, refine.tol=1e-6, nResample=100, bw=0.3, smooth=NULL, raf=c("NED","GKL","PWD","HD","SCHI2"), tau=1, subdivisions=1000, lambda.step=TRUE, sigma.step=TRUE, step=1, minw=0.04, nexp=1000, reparam=NULL, bootstrap=FALSE, bootstrap.lambda=NULL) {

#tuning.rho, tuning.chi # c1, c2
#lower,upper     # optimization interval
#n               # optimization grid
#max.it          # maxit in refinement
#refine.tol      # relative tolerance in refinement 
#nResample       # number of subsamples
## Parameters for Fully Itereated Weighted Likelihood and for One Step Weighted Likelihood 
#bw              # bandwidth
#smooth          # smooth parameter, if not NULL then bw = smooth*\hat{sd}
#raf             # residual adjustment function c("SCHI2", "NED", "HD", "GKL")
#tau             # parameter for the GKL and PD family
#subdivisions
#lambda.step
#sigma.step
#step            #step length in the oneWL function

#Used only in the oneWL function 
#minw            #below this value the weights are set to zero.
#nexp            #number of quantile points used to approximate the Expected Jacobian.
#reparam         #list for reparametrization of the sigma parameter.
#bootstrap       #to bootstrap the oneWL.
#bootstrap.lambda#to bootstrap the oneWL.  
  raf <- match.arg(raf)

  if (lower >= upper)
    stop("'lower' must be smaller than 'upper'")
  
  if (length(step) != 1 & length(step) != 3)
    stop("'step' must be a scalar or a vector of length 3")

  if (length(n) > 1 || n < 0)
    stop("'n' must be a positive scalar")
  n <- round(n)
  n <- ifelse(n < 1, 1, n)

  if (length(max.it) > 1 || max.it < 0)
    stop("'max.it' must be a positive scalar")
  max.it <- round(max.it)
  max.it <- ifelse(max.it < 1, 1, max.it)

  if (length(refine.tol) > 1 || refine.tol < 0)
    stop("'refine.tol' must be a positive scalar")

  if (length(bw) > 1 || bw < 0)
    stop("'bw' must be a positive scalar")
  if (!is.null(smooth) && (length(bw) > 1 || bw < 0))
    stop("'smooth' must be NULL or a positive scalar")
  
  if (length(subdivisions) > 1 || subdivisions < 0)
    stop("'subdivisions' must be a positive scalar")
  subdivisions <- round(subdivisions)
  subdivisions <- ifelse(subdivisions < 1, 1, subdivisions)

  if (!is.logical(lambda.step))
    stop("'lambda.step' must be a logical")

  if (!is.logical(sigma.step))
    stop("'sigma.step' must be a logical")

  if (length(minw) > 1 || minw > 1 | minw < 0)
    stop("'minw' must be a scalar in [0,1]")

  if (length(nexp) > 1 || nexp < 0)
    stop("'nexp' must be a positive scalar")
  nexp <- round(nexp)
  nexp <- ifelse(nexp < 1, 1, nexp)
  
############ EXAMPLE OF REPARAM ################################
#reparam.loggamma <- list(gam=function(sigma) sqrt(sigma),
#                         gaminv=function(gam) gam^2,
#                         delta=function(sigma) 2*sqrt(sigma))

  if (!is.null(reparam)) {
    badreparam <- !is.list(reparam) || !is.function(reparam$gam) || !is.function(reparam$gaminv) || !is.function(reparam$delta)
    if (badreparam)
      stop("'reparam' must be a list of 3 functions named 'gam', 'gaminv' and 'delta'")
  }

  if (!is.logical(bootstrap))
    stop("'bootstrap' must be a logical")

  res <- list(method=method, tuning.rho=tuning.rho, tuning.psi=tuning.psi, lower=lower, upper=upper, n=n, max.it=max.it, refine.tol=refine.tol, nResample=nResample, bw=bw, smooth=smooth, raf=raf, tau=tau, subdivisions=subdivisions, lambda.step=lambda.step, sigma.step=sigma.step,step=step,minw=minw, nexp=nexp, reparam=reparam, bootstrap=bootstrap, bootstrap.lambda=bootstrap.lambda)
  return(res)
}

#############################################################
#	loggammarob.QTau function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.QTau <- function(x, w=rep(1, length(x)), control=loggammarob.control()) {
  lgrid <- seq(control$lower, control$upper, length.out=control$n)
  if (!control$bootstrap) {
    resQ <- WQtau(ri=x,w=w,lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)
  } else {
    resQ <- WQtauboot(ri=x,w=w,lambda=control$bootstrap.lambda,step=0.25,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)
  }
  result <- list(mu=resQ$mu2, sigma=resQ$sig2, lambda=resQ$lam2, weights=w, iterations=NULL, error=NULL)
  return(result)
}

#############################################################
#	loggammarob.WQTau function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.WQTau <- function(x, w=rep(1, length(x)), control=loggammarob.control()) {
  lgrid <- seq(control$lower, control$upper, length.out=control$n)
  resQ <- WQtau(ri=x,w=w,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)  
  pp <- ppoints(length(x))
  ql <- qloggamma(p=pp, lambda=resQ$lam1)
  dl <- dloggamma(x=ql, lambda=resQ$lam1)
  vl <- pp*(1-pp)/dl^2
  wl <- 1/sqrt(vl)
  
  resWQ <- WQtau(ri=x,w=wl,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)

  result <- list(mu=resWQ$mu2, sigma=resWQ$sig2, lambda=resWQ$lam2, weights=wl, iterations=NULL, error=NULL, QTau=list(mu=resQ$mu2, sigma=resQ$sig2, lambda=resQ$lam2, weights=w, iterations=NULL, error=NULL))
  return(result)
}

#############################################################
#	loggammarob.WL function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.WL <- function(x, start=NULL, w=rep(1, length(x)), control=loggammarob.control()) {
  if (is.null(start)) {
    lgrid <- seq(control$lower, control$upper, length.out=control$n)
    resQ <- WQtau(ri=x,w=w,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)
    pp <- ppoints(length(x))
    ql <- qloggamma(p=pp, lambda=resQ$lam1)
    dl <- dloggamma(x=ql, lambda=resQ$lam1)
    vl <- pp*(1-pp)/dl^2
    wl <- 1/sqrt(vl)

    resWQ <- WQtau(ri=x,w=wl,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)
    start <- c(resWQ$mu2, resWQ$sig2, resWQ$lam2)
    QTau <- list(mu=resQ$mu2, sigma=resQ$sig2, lambda=resQ$lam2, weights=w, iterations=NULL, error=NULL)
    WQTau <- list(mu=resWQ$mu2, sigma=resWQ$sig2, lambda=resWQ$lam2, weights=wl, iterations=NULL, error=NULL)
  } else {
    if (!is.numeric(start))
      stop("'start' must be a numeric vector")
    if (length(start)!=3)
      stop("'start' must be a vector of length 3: mu, sigma2, lambda")
    QTau <- NULL
    WQTau <- NULL
  }
  
  if (!is.null(control$smooth))
    control$bw <- control$smooth*sqrt(start[2]) 
  resWL <- Disparity.WML.loggamma(y=x,mu0=start[1],sig0=start[2],lam0=start[3],lam.low=control$lower,lam.sup=control$upper,tol=control$refine.tol,maxit=control$max.it,lstep=control$lambda.step,sigstep=control$sigma.step,bw=control$bw,raf=control$raf,tau=control$tau,nmod=control$subdivisions,minw=control$minw)
  
  result <- list(mu=resWL$mu, sigma=resWL$sig, lambda=resWL$lam, weights=resWL$weights, iterations=resWL$nit, error=NULL, QTau=QTau, WQTau=WQTau)
  return(result)
}

#############################################################
#	loggammarob.ML function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.ML <- function(x, start=NULL, w=rep(1, length(x)), control=loggammarob.control()) {
  if (is.null(start)) {
    lgrid <- seq(control$lower, control$upper, length.out=control$n)
    resQ <- WQtau(ri=x,w=w,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)
    start <- c(resQ$mu2, resQ$sig2, resQ$lam2)    
  } else {
    if (!is.numeric(start))
      stop("'start' must be a numeric vector")
    if (length(start)!=3)
      stop("'start' must be a vector of length 3: mu, sigma2, lambda") 
  }
  resML <- WML.loggamma(y=x,wi=w,mu0=start[1],sig0=start[2],lam0=start[3],lam.low=control$lower,lam.sup=control$upper,tol=control$refine.tol,maxit=control$max.it,lstep=control$lambda.step,sigstep=control$sigma.step)
  
  result <- list(mu=resML$mu, sigma=resML$sig, lambda=resML$lam, weights=w, iterations=resML$nit, error=NULL)
  return(result)
}

#############################################################
#	loggammarob.oneWL function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.oneWL <- function(x, start=NULL, w=rep(1, length(x)), control=loggammarob.control()) {
  if (is.null(start)) {
    lgrid <- seq(control$lower, control$upper, length.out=control$n)
    if (!control$bootstrap) {
      resQ <- WQtau(ri=x,w=w,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)
      pp <- ppoints(length(x))
      ql <- qloggamma(p=pp, lambda=resQ$lam1)
      dl <- dloggamma(x=ql, lambda=resQ$lam1)
      vl <- pp*(1-pp)/dl^2
      wl <- 1/sqrt(vl)
    } else {
      resQ <- resWQ <- WQtauboot(ri=x,w=w,lambda=control$bootstrap.lambda,step=0.25,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)
      wl <- w
###    cat(resQ$mu1, resQ$sig1, resQ$lam1, '\n')
###    resQ <- list(); resQ$lam1 <- control$bootstrap.lambda
    }

    if (!control$bootstrap)
      resWQ <- WQtau(ri=x,w=wl,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol) 
  ## else
  ##   resWQ <- WQtauboot(ri=x,w=wl,lambda=control$bootstrap.lambda,step=0.25,lgrid=lgrid,c1=control$tuning.rho,c2=control$tuning.psi,N=control$nResample,maxit=control$max.it,tolr=control$refine.tol)

####    cat(resWQ$mu2, resWQ$sig2, resWQ$lam2, '\n')

    start <- c(resWQ$mu2, resWQ$sig2, resWQ$lam2)
    QTau <- list(mu=resQ$mu2, sigma=resQ$sig2, lambda=resQ$lam2, weights=w, iterations=NULL, error=NULL)
    WQTau <- list(mu=resWQ$mu2, sigma=resWQ$sig2, lambda=resWQ$lam2, weights=wl, iterations=NULL, error=NULL)    
  } else {
    if (!is.numeric(start))
      stop("'start' must be a numeric vector")
    if (length(start)!=3)
      stop("'start' must be a vector of length 3: mu, sigma2, lambda")    
    QTau <- NULL
    WQTau <- NULL
  }

  if (!is.null(control$smooth))
    control$bw <- control$smooth*sqrt(start[2]) 
  
  if (is.null(control$reparam))
    resOWL <- WMLone.loggamma(yi.sorted=x,mu0=start[1],sig0=start[2],lam0=start[3],bw=control$bw,raf=control$raf,tau=control$tau,nmod=control$subdivisions,step=control$step,minw=control$minw,nexp=control$nexp)
  else 
    resOWL <- WMLone.reparam.loggamma(yi.sorted=x,mu0=start[1],sig0=start[2],lam0=start[3],bw=control$bw,raf=control$raf,tau=control$tau,nmod=control$subdivisions,step=control$step,minw=control$minw,nexp=control$nexp,reparam=control$reparam)  
  result <- list(mu=resOWL$mu, sigma=resOWL$sig, lambda=resOWL$lam, weights=resOWL$weights, iterations=1, error=resOWL$error, step=control$step, QTau=QTau, WQTau=WQTau)
  return(result)
}

#############################################################
#	loggammarob.test function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

loggammarob.test <- function(x, mu=NULL, sigma=NULL, lambda=NULL, eta=NULL, type="Wald", conf.level = 0.95, prob=0.00001) {
  ## x an object from loggammarob
  if (!(x$method %in% c("WL", "oneWL", "ML")))
    stop("Inference is available only for 'WL', 'oneWL' and 'ML'")
  if (is.null(mu) & is.null(sigma) & is.null(lambda) & is.null(eta)) {
    mu <- 0
    sigma <- 1
    lambda <- 0
  }
  alternative <- "two.sided"  
  dname <- names(x$data)  
  type <- match.arg(type)
  components <- c(!is.null(mu), !is.null(sigma), !is.null(lambda))
  if (is.null(eta)) {
    estimate <- c(x$mu, x$sigma, x$lambda)[components]
    null.value <- c(mu, sigma, lambda)
    if (type=="Wald") {
      param <- estimate - null.value
      df <- length(param)      
      acov <- ACOV.ML.loggamma(sigma=x$sigma,lambda=x$lambda, prob=prob)
      if (df==1) {
        hatJ <- 1/diag(acov$cov)[components]
      } else {
        hatJ <- acov$Fisher.Info
        hatJ <- matrix(hatJ[as.logical(components%*%t(components))], nrow=length(param), byrow=TRUE)
      }
      wstat <- drop(sum(x$weights)*t(param)%*%hatJ%*%param)
      pval <- pchisq(q=wstat, df=df, ncp=0, lower.tail = FALSE, log.p = FALSE)
      estcov <- diag(acov$cov)[components]
      if (df==1)
        cint <- estimate+c(-1,1)*qnorm(1-(1-conf.level)/2)*sqrt(estcov)/sqrt(sum(x$weights))
      else       
        cint <- c(NA,NA)
    }
    names(wstat) <- "ww"
    names(df) <- "df"
    names(null.value) <- c("mean", "scale", "shape")[components]    
  } else {
    if (any(components))
      warning("test with 'eta' together with other parameters is not implemented, only the test on 'eta' is performed")
    eps <- 0.0001
    estimate <- x$eta
    null.value <- eta
    if (type=="Wald") {
      param <- estimate - null.value
      if (x$lambda < eps)
        warning("test and confidence intervals for 'eta' is not implemented for negative lambda")
      hatvar <- AVAR.ML.eta.loggamma(mu=x$mu, sigma=x$sigma, lambda=x$lambda, prob=prob, eps=eps)
      wstat <- sum(x$weights)*param^2/hatvar
      df <- length(param)
      pval <- pchisq(q=wstat, df=df, ncp=0, lower.tail = FALSE, log.p = FALSE)
      cint <- estimate+c(-1,1)*qnorm(1-(1-conf.level)/2)*sqrt(hatvar)/sqrt(sum(x$weights))
    }
    names(wstat) <- "ww"
    names(df) <- "df"
    names(null.value) <- "mean(exp(X))"        
  }
  
  type <- paste(ifelse(x$method=="ML","Classic","Weighted"), type, "Test based on", x$method, sep=" ")
  attr(cint,"conf.level") <- conf.level
  rval <- list(statistic=wstat, parameter=df, p.value=pval,
    conf.int=cint, estimate=estimate, null.value=null.value,
    alternative=alternative,
    method=type, data.name=dname)
    class(rval) <- "htest"
  return(rval)
}

#############################################################
#	print.loggammarob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

######## PRINT
print.loggammarob <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
  cat("Location: ", format(x$mu, digits=digits))
  cat("  Scale: ", format(x$sigma, digits=digits))
  cat("  Shape: ", format(x$lambda, digits=digits))
  cat("  E(exp(X)): ", format(x$eta, digits=digits))  
  cat("\n")
  invisible(x)
}

#############################################################
#	summary.loggammarob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

########### SUMMARY
summary.loggammarob <- function(object, p=NULL, conf.level=0.95, prob=0.00001, ...) {
  if ((object$method %in% c("WL", "oneWL", "ML"))) {
    p <- p[ifelse(p > 0 & p < 1, TRUE, FALSE)]
    nrep <- 4+length(p)
    if (length(conf.level) == 1)
      conf.level <- rep(conf.level, nrep)  
    if (length(conf.level) != nrep)
      stop("'conf.level' must have length equal to 1 or 4+length(p)")

    z <- qnorm(1-(1-conf.level)/2)
    if (is.finite(object$eta)) {
      asd <- sqrt(diag(ACOV.ML.loggamma(sigma=object$sigma,lambda=object$lambda,prob=prob)$cov))
      asdeta <- sqrt(AVAR.ML.eta.loggamma(mu=object$mu, sigma=object$sigma, lambda=object$lambda, prob=prob, eps=0.0001, npoints=100000))
    } else {
      asd <- rep(NA,3)
      asdeta <- NA
    }
    object$muconf.int <- object$mu+c(-1,1)*z[1]*asd[1]/sqrt(sum(object$weights))
    object$sigmaconf.int <- object$sigma+c(-1,1)*z[2]*asd[2]/sqrt(sum(object$weights))
    object$lambdaconf.int <- object$lambda+c(-1,1)*z[3]*asd[3]/sqrt(sum(object$weights))
    object$etaconf.int <- object$eta+c(-1,1)*z[4]*asdeta/sqrt(sum(object$weights))

    attr(object$muconf.int, "conf.level") <- conf.level[1]
    attr(object$sigmaconf.int, "conf.level") <- conf.level[2]
    attr(object$lambdaconf.int, "conf.level") <- conf.level[3]
    attr(object$etaconf.int, "conf.level") <- conf.level[4]

    object$muse <- asd[1]/sqrt(sum(object$weights))
    object$sigmase <- asd[2]/sqrt(sum(object$weights))
    object$lambdase <- asd[3]/sqrt(sum(object$weights))
    object$etase <- asdeta/sqrt(sum(object$weights))

    object$p <- p
    
    if (!is.null(p)) {
      object$q <- asdq <- rep(NA, length(p))
      object$qconf.int <- matrix(NA, nrow=length(p), ncol=2)
      for (i in 1:length(p)) {
        object$q[i] <- qloggamma(p=p[i], mu=object$mu, sigma=object$sigma, lambda=object$lambda)
        if (is.finite(object$eta))
          asdq[i] <- sqrt(AVAR.ML.quantile.loggamma(p=p[i], mu=object$mu, sigma=object$sigma, lambda=object$lambda, prob=prob))
        else
          asdq[i] <- NA
        object$qconf.int[i,] <- object$q[i]+c(-1,1)*z[4+i]*asdq[i]/sqrt(sum(object$weights))
      }
#####      object$qconf.int <- drop(qconf.int)
      attr(object$qconf.int, "conf.level") <- conf.level[4+(1:length(p))]
      object$qse <- asdq/sqrt(sum(object$weights))
    }
  }

  object$call <- match.call()
  class(object) <- "summary.loggammarob"
  return(object)
}

#############################################################
#	print.summary.loggammarob function
#	Author: C. Agostinelli, A. Marazzi,
#               V.J. Yohai and A. Randriamiharisoa
#	Maintainer e-mail: claudio@unive.it
#	Date: January, 1, 2013
#	Version: 0.1
#	Copyright (C) 2013 C. Agostinelli A. Marazzi,
#                  V.J. Yohai and A. Randriamiharisoa
#############################################################

#PRINT.SUMMARY
print.summary.loggammarob <- function(x, digits = max(3, getOption("digits") - 3), ...) {
  cat("\nCall:\n",deparse(x$call),"\n\n",sep="")
  cat("Location: ", format(x$mu, digits=digits))
  if ((x$method %in% c("WL", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$muse, digits=digits), "\n")
    cat("(", format(x$muconf.int[1L], digits=digits), ", ",
        format(x$muconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$muconf.int, "conf.level")),
        "percent confidence interval\n\n")
  } else {
    cat("\n")
  }
  cat("Scale: ", format(x$sigma, digits=digits))
  if ((x$method %in% c("WL", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$sigmase, digits=digits), "\n")    
    cat("(", format(x$sigmaconf.int[1L], digits=digits), ", ",
        format(x$sigmaconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$sigmaconf.int, "conf.level")),
        "percent confidence interval\n\n")
  } else {
    cat("\n")
  }
  cat("Shape: ", format(x$lambda, digits=digits))
  if ((x$method %in% c("WL", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$lambdase, digits=digits), "\n")    
    cat("(", format(x$lambdaconf.int[1L], digits=digits), ", ",
        format(x$lambdaconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$lambdaconf.int, "conf.level")),
        "percent confidence interval\n\n")
  } else {
    cat("\n")
  }
  cat("Mean(exp(X)): ", format(x$eta, digits=digits))
  if ((x$method %in% c("WL", "oneWL", "ML"))) {
    cat(" s.e. ", format(x$etase, digits=digits), "\n")    
    cat("(", format(x$etaconf.int[1L], digits=digits), ", ",
        format(x$etaconf.int[2L], digits=digits), ") \n")
    cat(format(100 * attr(x$etaconf.int, "conf.level")),
        "percent confidence interval\n\n\n")
  } else {
    cat("\n")
  }

######### Quantiles estimation and confidence interval
  if (!is.null(x$p)) {
    for (i in 1:length(x$p)) {
      cat("Quantile of order ", format(x$p[i], digits=digits), ": ", format(x$q[i], digits=digits))
      if ((x$method %in% c("WL", "oneWL", "ML"))) {
        cat(" s.e. ", format(x$qse[i], digits=digits), "\n")    
        cat("(", format(x$qconf.int[i,1L], digits=digits), ", ",
        format(x$qconf.int[i,2L], digits=digits), ") \n")
        cat(format(100 * attr(x$qconf.int, "conf.level")[i]),
          "percent confidence interval\n\n\n")
      } else {
        cat("\n")
      }
    }
  }

  if ((x$method %in% c("WL", "oneWL", "ML")))
    cat("\n\n")
  else
    cat("\n\nConfidence intervals and Standard Errors are available only for 'WL', 'oneWL' and 'ML'\n\n\n") 
  
  if (x$method %in% c("WL", "oneWL", "WQTau"))
    summarizeRobWeights(x$weights, digits = digits, ...)
  invisible(x)
}
