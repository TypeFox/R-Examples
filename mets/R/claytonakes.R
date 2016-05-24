##' Clayton-Oakes frailty model
##'
##' @title Clayton-Oakes model with piece-wise constant hazards
##' @param formula formula specifying the marginal proportional (piecewise constant) hazard structure with the right-hand-side being a survival object (Surv) specifying the entry time (optional), the follow-up time, and event/censoring status at follow-up. The clustering can be specified using the special function \code{cluster} (see example below).
##' @param data Data frame
##' @param cluster Variable defining the clustering (if not given in the formula)
##' @param var.formula Formula specifying the variance component structure (if not given via the cluster special function in the formula) using a linear model with log-link.
##' @param cuts Cut points defining the piecewise constant hazard
##' @param type when equal to \code{two.stage}, the Clayton-Oakes-Glidden estimator will be calculated via the \code{timereg} package
##' @param start Optional starting values
##' @param control Control parameters to the optimization routine
##' @param var.invlink Inverse link function for variance structure model
##' @param ... Additional arguments
##' @author Klaus K. Holst
##' @examples
##' set.seed(1)
##' d <- subset(simClaytonOakes(500,4,2,1,stoptime=2,left=2),!truncated)
##' e <- ClaytonOakes(Surv(lefttime,time,status)~x1+cluster(~1,cluster),
##'                   cuts=c(0,0.5,1,2),data=d)
##' e
##' 
##' 
##' d2 <- simClaytonOakes(500,4,2,1,stoptime=2,left=0)
##' d2$z <- rep(1,nrow(d2)); d2$z[d2$cluster%in%sample(d2$cluster,100)] <- 0
##' ## Marginal=Cox Proportional Hazards model:
##' ts <- ClaytonOakes(Surv(time,status)~prop(x1)+cluster(~1,cluster),
##'                    data=d2,type="two.stage")
##' ## Marginal=Aalens additive model:
##' ts2 <- ClaytonOakes(Surv(time,status)~x1+cluster(~1,cluster),
##'                     data=d2,type="two.stage")
##' ## Marginal=Piecewise constant:
##' e2 <- ClaytonOakes(Surv(time,status)~x1+cluster(~-1+factor(z),cluster),
##'                    cuts=c(0,0.5,1,2),data=d2)
##' e2
##' plot(ts)
##' plot(e2,add=TRUE)
##' 
##' e3 <- ClaytonOakes(Surv(time,status)~x1+cluster(~1,cluster),cuts=c(0,0.5,1,2),
##'                    data=d,var.invlink=identity)
##' e3
##' @export
ClaytonOakes <- function(formula,data=parent.frame(),cluster,var.formula=~1,cuts=NULL,type="piecewise",start,control=list(),var.invlink=exp,...) {

  mycall <- match.call()
  dots <- list(...)
  formulaId <- Specials(formula,"cluster") 
  formulaStrata <- Specials(formula,"strata")
  formulaSt <- "~."
  formulaProp <- Specials(formula,"prop")
  if (!is.null(formulaId)) {
    var.formulaId <- ~1
    if (length(formulaId)>1) {
      var.formula <- as.formula(formulaId[[1]])
      formulaId <- formulaId[[2]]
    }
    cluster <- formulaId
    mycall$cluster <- cluster
    formulaSt <- paste(formulaSt,paste("-cluster(",paste(var.formula,collapse=""),
                                       ",",formulaId,")"))
  }
  formulaSt <- paste(formulaSt,paste("-strata(",paste(formulaStrata,collapse="+"),")"))
  formula <- update(formula,formulaSt)
  
  if (!is.null(formulaStrata)) {
    strata <- formulaStrata
    mycall$strata <- strata
  }
  if (missing(cluster)) stop("Missing 'cluster' variable")
  ngamma <- 0
  data <- data[order(data[,cluster]),]
  Z <- model.matrix(var.formula,data)
  ngamma <- ncol(Z)

  if (type!="piecewise") {
    timeregmod <- ifelse(length(formulaProp)>0,"cox.aalen","aalen")
    if (is.null(dots$robust)) dots$robust <- 0
    args <- c(list(formula=formula,data=data,max.clust=NULL,clusters=data[,cluster]),dots)
    marg <- do.call(timeregmod, args)    
    return(two.stage(marg,data=data,theta.des=Z,var.link=1,...))
  }
  
  
  timevar <- terms(formula)[[2]]
  if (is.call(timevar)) {
    delayedentry <- (length(timevar)==4)*1
    entry <- NULL
    if (delayedentry==1)
      entry <- as.character(timevar[[2]])
    causes <- timevar[[3+delayedentry]]
    timevar <- timevar[[2+delayedentry]]
  }  
  timevar <- as.character(timevar)
  causes <- as.character(causes)
  covars <- as.character(attributes(terms(formula))$variables)[-(1:2)]
  X <- NULL
  nbeta <- 0  
  if (length(covars)>0) {
##    X <- model.matrix(as.formula(paste("~-1+",paste(covars,collapse="+"))),data)
    X <- model.matrix(update(formula,.~.+1),data)[,-1,drop=FALSE]
    nbeta <- ncol(X)
  }
  
  if (is.data.frame(data)) {
    mydata <- data.frame(T=data[,timevar],status=data[,causes],cluster=data[,cluster],entry=0)
    if (!is.null(entry)) {
      mydata$entry <- data[,entry]
    }
  } else {
    mydata <- data.frame(T=get(timevar,envir=data),status=get(causes,envir=data),cluster=get(cluster,envir=data),entry=0)
    if (!is.null(entry))      
      mydata$entry <- get(entry,envir=data)
  }
  if (is.null(cuts)) {
    cuts <- c(0,max(mydata$T))
  }
  if (max(mydata$T)>tail(cuts,1)) stop("Interval does not embed time observations")
  if (any(with(mydata, T<entry))) stop("Entry time occuring after event")

  ucluster <- unique(mydata$cluster)
  
  npar <- length(cuts)-1
  if (!is.null(X)) npar <- npar+ncol(X)
  npar <- npar+ncol(Z)
  p0 <- rep(0.1,npar)
  if (!missing(start)) p0 <- c(start,rep(0,max(0,length(npar)-length(start))))
  if (!is.null(dots$theta)) p0[1] <- dots$theta

  invlinkname <- as.character(substitute(var.invlink))

  
  obj <- function(p) {
    varpar <- p[seq(ngamma)]
    p <- p[-seq(ngamma)]
    ##    theta0 <- rep(exp(varpar),length(ucluster));
    theta0 <- Z%*%varpar
    if (exists(invlinkname)) theta0 <- var.invlink(theta0)
    multhaz <- rep(1,nrow(mydata))
    if (!is.null(X)) {
      nbeta <- ncol(X)
      beta <- p[seq(nbeta)]
      p <- p[-seq(nbeta)]
      multhaz <- exp(X%*%beta)
    }
    if (!is.null(dots$theta)) theta0 <- cbind(rep(dots$theta,length(theta0)))
    res <- .Call("claytonoakes",
           ds=mydata$status,ts=mydata$T,es=mydata$entry,
           allcs=mydata$cluster,cs=ucluster, cuts=cuts,
                 hs=exp(p),mulths=multhaz,
                 var=theta0)$logLik    
    return(-res)
  }
  opt <- tryCatch(nlminb(p0,obj,control=control),error=function(x) NULL)
  if (is.null(opt)) stop("Critical optmization problem")
  if (any(is.na(opt$par)) | any(!is.finite(opt$par)) | any(is.nan(opt$par)) ) {
    V <- matrix(NA,length(p0),length(p0))
  } else {    
    I <- numDeriv::hessian(obj,opt$par)
    ee <- tryCatch(eigen(I),error=function(x) NULL); 
    if (!is.null(ee)) {
      threshold <- 1e-12
      idx <- ee$values>threshold
      ee$values[idx] <- 1/ee$values[idx];
      if (!all(idx))
        ee$values[!idx] <- 0
      V <- with(ee, vectors%*%diag(values)%*%t(vectors))
    } else {
      V <- matrix(NA,length(p0),length(p0))
    }
  }
  res <- list(coef=opt$par,vcov=V,cuts=cuts,nbeta=nbeta,ngamma=ngamma,betanames=colnames(X),gammanames=colnames(Z),opt=opt,invlink=var.invlink,invlinkname=invlinkname)
  class(res) <- "claytonoakes"
  return(res)
}

##################################################

##' @export
print.claytonoakes <- function(x,...) {
  print(summary(x))
}

##' @export
print.summary.claytonoakes <- function(x,...) {
  printCoefmat(x$coef[,c(1,3,4)],...)
  cat("\nDependence parameters:\n")
  printCoefmat(x$var,...)
  invisible(x)
}

##' @export
summary.claytonoakes <- function(object,...) {
  mycoef <- matrix(nrow=length(object$coef),ncol=4)
  mycoef[,1:2] <- cbind(object$coef,sqrt(diag(object$vcov)))
  mycoef[,3:4] <- cbind(mycoef[,1]-qnorm(0.975)*mycoef[,2],mycoef[,1]+qnorm(0.975)*mycoef[,2])
  colnames(mycoef) <- c("Estimate","Std.Err","2.5%","97.5%")
  if (length(object$cuts))
  cutnames <- levels(cut(0,breaks=object$cuts))
  varname <- switch(object$invlinkname,exp="log-Var:",identity="Var:",paste("inv",object$invlinkname,"-Var:",sep=""))
  rownames(mycoef) <- c(paste(varname,object$gammanames,sep=""),object$betanames,cutnames)
  mycoef[-seq(object$ngamma),] <- exp(mycoef[-seq(object$ngamma),])
  varcoef <- object$invlink(mycoef[seq(object$ngamma),c(1,3,4),drop=FALSE])
  rownames(varcoef) <- object$gammanames
  varcoef <- cbind(varcoef,1/(1+2/varcoef))
  colnames(varcoef)[c(1,4)] <- c("Variance","Kendall's tau")  
  res <- list(coef=mycoef,var=varcoef)
  class(res) <- "summary.claytonoakes"
  res
}

##' @export
plot.claytonoakes <- function(x,chaz=TRUE,add=!is.null(dev.list()),col="darkblue",...) {
  haz <- summary(x)$coef[-seq(x$nbeta+x$ngamma),,drop=FALSE]
  t <- x$cuts
  L <- approxfun(t,f=1,cumsum(c(0,haz[,1]*diff(t))),method="linear")
  if (add) {
    lines(t,L(t),col=col,...)
  } else {
    plot(t,L(t),type="l",col=col,...)
  }
  invisible(x)  
}

predict.claytonoakes <- function(x,...) {

}
