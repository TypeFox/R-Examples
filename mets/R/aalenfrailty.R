##' Additive hazards model with (gamma) frailty
##'
##' Aalen frailty model
##' @title Aalen frailty model
##' @param time Time variable
##' @param status Status variable (0,1)
##' @param X Covariate design matrix
##' @param id cluster variable
##' @param theta list of thetas (returns score evaluated here), or
##' starting point for optimization (defaults to magic number 0.1)
##' @param B (optional) Cumulative coefficients (update theta by fixing B)
##' @param ... Additional arguments to lower level functions
##' @return Parameter estimates
##' @author Klaus K. Holst
##' @export
##' @examples
##' dd <- simAalenFrailty(5000)
##' f <- ~1##+x
##' X <- model.matrix(f,dd) ## design matrix for non-parametric terms
##' system.time(out<-aalen(update(f,Surv(time,status)~.),dd,n.sim=0,robust=0))
##' dix <- which(dd$status==1)
##' t1 <- system.time(bb <- .Call("Bhat",as.integer(dd$status),
##'                               X,0.2,as.integer(dd$id),NULL,NULL,-1,
##'                               package="mets"))
##' spec <- 1
##' ##plot(out,spec=spec)
##' ## plot(dd$time[dix],bb$B2[,spec],col="red",type="s",
##' ##      ylim=c(0,max(dd$time)*c(beta0,beta)[spec]))
##' ## abline(a=0,b=c(beta0,beta)[spec])
##' ##'
##' 
##' \dontrun{
##' thetas <- seq(0.1,2,length.out=10)
##' Us <- unlist(aalenfrailty(dd$time,dd$status,X,dd$id,as.list(thetas)))
##' ##plot(thetas,Us,type="l",ylim=c(-.5,1)); abline(h=0,lty=2); abline(v=theta,lty=2)
##' op <- aalenfrailty(dd$time,dd$status,X,dd$id)
##' op
##' }
aalenfrailty <- function(time,status,X,id,theta,B=NULL,...) {  
  dix <- which(status==1)
  cc <- cluster.index(id)
  ncluster <- length(cc$clusters)
  U <- function(theta,indiv=FALSE,Bhat=FALSE) {
    if (is.null(B)) {
        BB <- .Call("Bhat",as.integer(status),X,theta,as.integer(cc$clusters),cc$idclust,as.integer(cc$cluster.size))$B2
    } else {
        BB <- B*time[dix]
    }
    Hij0 <- apply(X[dix,,drop=FALSE]*BB,1,sum)
    Hij <- Cpred(cbind(time[dix],Hij0),time)[,2,drop=FALSE]
##    if (is.na(Hij[1])) browser()
    res <- .Call("Uhat",as.integer(status),Hij,theta,cc$idclust,as.integer(cc$cluster.size))
    if (!indiv) res <- mean(res,na.rm=TRUE)
    if (Bhat) attributes(res)$B <- BB
    return(res)
  }
  if (missing(theta)) theta <- 0.1
  if (is.list(theta)) {
      cc <- lapply(theta,function(x) U(x,Bhat=TRUE,...))
      BB <- Reduce("cbind",lapply(cc,function(x) attributes(x)$B))
      UU <- unlist(lapply(cc,function(x) x[1]))
      res <- list(U=UU, B=BB, time=time, dix=dix, X=X, id=id, status=status)
      return(res)
  }
  op <- nlminb(theta,function(x) U(x)^2)
  uu <- U(op$par,TRUE)
  du <- numDeriv::grad(U,op$par)
  return(list(theta=op$par, sd=(mean(uu^2)/du^2/ncluster)^0.5))
}
                  

##' Simulate observations from Aalen Frailty model with Gamma
##' distributed frailty and constant intensity.
##'
##' @title Simulate from the Aalen Frailty model
##' @param n Number of observations in each cluster
##' @param theta Dependence paramter (variance of frailty)
##' @param K Number of clusters
##' @param beta0 Baseline (intercept)
##' @param beta Effect (log hazard ratio) of covariate
##' @param cens Censoring rate
##' @param cuts time cuts
##' @param ... Additional arguments
##' @author Klaus K. Holst
##' @export
simAalenFrailty <- function(n=5e3,theta=0.3,K=2,beta0=1.5,beta=1,cens=1.5,cuts=0,...) {
    ## beta0 (constant baseline intensity)
    ## beta (covariate effect)
    if (length(beta0)!=length(cuts)) 
        stop("Number of time-intervals (cuts) does not agree with number of rate parameters (beta0)")  
    cuts <- c(cuts,Inf)
    id <- rep(seq(n/K),each=K) ## Cluster indicator
    x <- rbinom(n,1,0.5) ## Binary covariate
    Z <- rep(rgamma(n/K,1/theta,1/theta),each=K) ## Frailty, mean 1, var theta
    Ai <- function() {
        vals <- matrix(0,ncol=length(beta0),nrow=n)
        ival <- numeric(n)
        for (i in seq(length(beta0))) {
            u <- -log(runif(n)) ##rexp(n,1)
            vals[,i] <-  cuts[i] + u/(beta0[i]+beta*x)/Z
            idx <- which(vals[,i]<=cuts[i+1] & ival==0)
            ival[idx] <- vals[idx,i]
        }
        ival
    }
    dat <- data.frame(time=Ai(), x=x, status=1, id=id, Z=Z)
    if (cens==0) cens <- Inf else cens <- -log(runif(n))/cens
    dat$status <- (dat$time<=cens)*1
    dat$time <- apply(cbind(dat$time,cens),1,min)
    dat <- dat[order(dat$time),] ## order after event/censoring time
    return(dat)
}
