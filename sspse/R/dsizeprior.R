#' Prior distributions for the size of a hidden population
#' 
#' \code{\link{dsizeprior}} computes the prior distribution of the population
#' size of a hidden population. The prior is intended to be used in Bayesian
#' inference for the population size based on data collected by Respondent
#' Driven Sampling, but can be used with any Bayesian method to estimate
#' population size.
#' 
#' 
#' @param n count; the sample size.
#' @param type character; the type of parametric distribution to use for the
#' prior on population size. The options are \code{"beta"} (for a Beta-type
#' prior on the sample proportion (i.e. \eqn{n/N}), \code{"nbinom"}
#' (Negative-Binomial), \code{"pln"} (Poisson-log-normal), \code{"flat"}
#' (uniform), \code{continuous} (the continuous version of the Beta-type prior
#' on the sample proportion). The last option is \code{"supplied"} which
#' enables a numeric prior to be specified. See the argument \code{supplied}
#' for the format of the information. The default \code{type} is \code{beta}.
#' @param mean.prior.size scalar; A hyperparameter being the mean of the prior
#' distribution on the population size.
#' @param sd.prior.size scalar; A hyperparameter being the standard deviation
#' of the prior distribution on the population size.
#' @param mode.prior.sample.proportion scalar; A hyperparameter being the mode
#' of the prior distribution on the sample proportion \eqn{n/N}.
#' @param median.prior.sample.proportion scalar; A hyperparameter being the
#' median of the prior distribution on the sample proportion \eqn{n/N}.
#' @param median.prior.size scalar; A hyperparameter being the mode of the
#' prior distribution on the population size.
#' @param mode.prior.size scalar; A hyperparameter being the mode of the prior
#' distribution on the population size.
#' @param quartiles.prior.size vector of length 2; A pair of hyperparameters
#' being the lower and upper quartiles of the prior distribution on the
#' population size. For example, \cr
#' \code{quartiles.prior.size=c(1000,4000)}
#' corresponds to a prior where the lower quartile (25\%) is 1000 and the upper
#' (75\%) is 4000.
#' @param effective.prior.df scalar; A hyperparameter being the effective
#' number of samples worth of information represented in the prior distribution
#' on the population size. By default this is 1, but it can be greater (or
#' less!) to allow for different levels of uncertainty.
#' @param alpha scalar; A hyperparameter being the first parameter of the Beta
#' prior model for the sample proportion. By default this is NULL, meaning that
#' 1 is chosen. it can be any value at least 1 to allow for different levels of
#' uncertainty.
#' @param beta scalar; A hyperparameter being the second parameter of the Beta
#' prior model for the sample proportion. By default this is NULL, meaning that
#' 1 is chosen. it can be any value at least 1 to allow for different levels of
#' uncertainty.
#' @param maxN integer; maximum possible population size. By default this is
#' determined from an upper quantile of the prior distribution.
#' @param log logical; return the prior or the the logarithm of the prior.
#' @param maxbeta integer; maximum beta in the prior for population size. By
#' default this is determined to ensure numerical stability.
#' @param maxNmax integer; maximum possible population size. By default this is
#' determined to ensure numerical stability.
#' @param supplied list; If the argument \code{type="supplied"} then this
#' should be a list object, typically of class \code{sspse}. It is primarily
#' used to pass the posterior sample from a separate \code{size} call for use
#' as the prior to this call. Essentially, it must have two components named
#' \code{maxN} and \code{sample}. \code{maxN} is the maximum population
#' envisaged and \code{sample} is random sample from the prior distribution.
#' @param verbose logical; if this is \code{TRUE}, the program will print out
#' additional information, including goodness of fit statistics.
#' @return \code{\link{dsizeprior}} returns a list consisting of the following
#' elements: \item{x}{vector; vector of degrees \code{1:N} at which the prior
#' PMF is computed.} \item{lpriorm}{vector; vector of probabilities
#' corresponding to the values in \code{x}.} \item{N}{scalar; a starting value
#' for the population size computed from the prior.} \item{maxN}{integer;
#' maximum possible population size. By default this is determined from an
#' upper quantile of the prior distribution.} \item{mean.prior.size}{scalar; A
#' hyperparameter being the mean of the prior distribution on the population
#' size.} 
#' \item{mode.prior.size}{scalar; A hyperparameter being the mode of the prior
#' distribution on the population size.} \item{effective.prior.df}{scalar; A
#' hyperparameter being the effective number of samples worth of information
#' represented in the prior distribution on the population size. By default
#' this is 1, but it can be greater (or less!) to allow for different levels of
#' uncertainty.} \item{mode.prior.sample.proportion}{scalar; A hyperparameter
#' being the mode of the prior distribution on the sample proportion
#' \eqn{n/N}.} \item{median.prior.size}{scalar; A hyperparameter being the mode
#' of the prior distribution on the population size.} \item{beta}{scalar; A
#' hyperparameter being the second parameter of the Beta distribution that is a
#' component of the prior distribution on the sample proportion \eqn{n/N}.}
#' \item{type}{character; the type of parametric distribution to use for the
#' prior on population size. The possible values are \code{beta} (for a Beta
#' prior on the sample proportion (i.e. \eqn{n/N}), \code{nbinom}
#' (Negative-Binomial), \code{pln} (Poisson-log-normal), \code{flat} (uniform),
#' and \code{continuous} (the continuous version of the Beta prior on the
#' sample proportion. The default is \code{beta}.}
#' @section Details on priors: The best way to specify the prior is via the
#' hyperparameter \code{mode.prior.size} which specifies the mode of the prior
#' distribution on the population size. You can alternatively specify the
#' hyperparameter \code{median.prior.size} which specifies the median of the
#' prior distribution on the population size, or \code{mode.prior.sample
#' proportion} which specifies the mode of the prior distribution on the
#' proportion of the population size in the sample.
#' @seealso network, statnet, degreenet
#' @references
#' 
#' Gile, Krista J. (2008) \emph{Inference from Partially-Observed Network
#' Data}, Ph.D. Thesis, Department of Statistics, University of Washington.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2010) \emph{Respondent-Driven
#' Sampling: An Assessment of Current Methodology}, Sociological Methodology
#' 40, 285-327.
#' 
#' Gile, Krista J. and Handcock, Mark S. (2014) \pkg{sspse}: Estimating Hidden 
#' Population Size using Respondent Driven Sampling Data
#' R package, Los Angeles, CA.  Version 0.5, \url{http://hpmrg.org}.
#' 
#' Handcock MS (2003).  \pkg{degreenet}: Models for Skewed Count Distributions
#' Relevant to Networks.  Statnet Project, Seattle, WA.  Version 1.2,
#' \url{http://statnetproject.org}.
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2014)
#' \emph{Estimating Hidden Population Size using Respondent-Driven Sampling
#' Data}, Electronic Journal of Statistics, 8, 1, 1491-1521
#' 
#' Handcock, Mark S., Gile, Krista J. and Mar, Corinne M. (2015)
#' \emph{Estimating the Size of Populations at High Risk for HIV using Respondent-Driven 
#' Sampling Data}, Biometrics.
#' @keywords models
#' @examples
#' 
#' dsizeprior(n=100,
#'            type="beta",
#'            mode.prior.size=1000)
#' 
#' @export dsizeprior
dsizeprior<-function(n,
		  type=c("beta","nbinom","pln","flat","continuous","supplied"),
		  mean.prior.size=NULL, sd.prior.size=NULL,
		  mode.prior.sample.proportion=NULL,
		  median.prior.sample.proportion=NULL,
		  median.prior.size=NULL,
		  mode.prior.size=NULL,
		  quartiles.prior.size=NULL,
		  effective.prior.df=1,
		  alpha=NULL,
		  beta=NULL,
                  maxN=NULL,
                  log=FALSE,
                  maxbeta=120,
                  maxNmax=200000,
                  supplied=list(maxN=maxN),
                  verbose=TRUE){
  priorsizedistribution=match.arg(type)
  N <- NULL
  maxN.set <- maxN
  lfn <- function(x,beta,n,effective.prior.df,alpha){
   a=effective.prior.df*(log(n)+lgamma(alpha+beta)-lgamma(alpha)-lgamma(beta)+(beta-1)*log(x-n) - (alpha+beta)*log(x) )
   mina <- min(a,na.rm=TRUE)
   a[is.na(a)] <- mina
   a[is.infinite(a)] <- mina
   a <- a - mina
   a
  }
  dfn <- function(alpha,beta,x,n,effective.prior.df){
   lpriorm <- lfn(x+0.5,beta,n,effective.prior.df,alpha)
   priorm <- exp(lpriorm)
   priorm/sum(priorm,na.rm=TRUE)
  }
  lpriorm <- switch(priorsizedistribution,
    nbinom={
      if(is.null(sd.prior.size)){sd.prior.size <- mean.prior.size}
      if(is.null(maxN)){
        maxN <- min(maxNmax,ceiling(qnbinom(p=0.995,
                    mu=mean.prior.size, prob=mean.prior.size/(sd.prior.size^2))))
      }
      if(is.null(N)){
        maxN <- min(maxNmax,ceiling(qnbinom(p=0.5,
                    mu=mean.prior.size, prob=mean.prior.size/(sd.prior.size^2))))
      }
      lpriorm <- dnbinom(x=n:maxN,
                           mu=mean.prior.size, prob=mean.prior.size/(sd.prior.size^2),
                           log=log)
      if(is.null(median.prior.size)) median.prior.size <- maxN/2
      lpriorm
     },
    flat={
      if(is.null(maxN)){
        maxN <- 10*n
      }
      if(is.null(N)){
        N <- 0.5*maxN
      }
#     if(is.null(mode.prior.size)){
#       mode.prior.size <- 0
#     }
      lpriorm <- rep(1/(maxN-n+1),maxN-n+1)
      if(log){
       lpriorm <- log(lpriorm)
#     }else{
#      if(mode.prior.size > n){
#       lpriorm[1:(mode.prior.size-n)] <- 0
#      }
      }
#     if(is.null(median.prior.size)) median.prior.size <- maxN/2
      lpriorm
     },
    continuous={
     if(!is.null(mode.prior.sample.proportion)){
      beta <- 2/mode.prior.sample.proportion - 1
     }
     if(!is.null(median.prior.sample.proportion)){
      beta <- -log(2)/log(1-median.prior.sample.proportion)
     }
     if(!is.null(median.prior.size)){
      beta <- -log(2)/log(1-n/median.prior.size)
     }
     if(!is.null(mean.prior.size)){
      beta <- mean.prior.size/n - 1
     }
     if(!is.null(mode.prior.size)){
      beta <- 2*mode.prior.size/n - 1
     }
     if(is.null(beta)){
       warning("No prior information about the population size was specified! Using a prior mode of twice the sample size. Please specify prior information!", call. = FALSE)
       beta <- 3
     }
     if(is.infinite(beta) | is.na(beta)){
       stop("The sample size is incompatible with the specified prior information about the population size.", call. = FALSE)
     }
     median.prior.size <- n/(1-0.5^(1/beta))
     mode.prior.size <- n*(beta+1)/2
     mode.prior.sample.proportion <- 2/(beta+1)
     if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
     if(is.null(N)){N <- min(maxNmax,ceiling( n/(1-0.5^(1/beta)) ))}
     x <- (n:maxN)
     lpriorm <- log(beta*n)+(beta-1)*log(x-n+1)-(beta+1)*log(x+1)
     if(!log){
      lpriorm <- exp(lpriorm)
     }
     lpriorm
     },
    beta={
     if(is.null(alpha) | is.null(beta)){
     if(!is.null(mode.prior.sample.proportion)){
      beta <- 2/mode.prior.sample.proportion - 1
     }
     if(!is.null(median.prior.sample.proportion)){
      beta <- -log(2)/log(1-median.prior.sample.proportion)
     }
     if(!is.null(median.prior.size)){
      if(median.prior.size < n){median.prior.size = n}
#     if(median.prior.size < 750){effective.prior.df=max(effective.prior.df,3)}
      beta <- -log(2)/log(1-n/median.prior.size)
      if(is.null(alpha)) alpha=median.prior.size/(median.prior.size-n)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn1 <- function(beta,x,n,median.prior.size,effective.prior.df,alpha){
       priorm <- dfn(alpha,beta,x,n,effective.prior.df)
       abs(median.prior.size - 
        (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)] ) 
      }
      if(is.null(maxN.set)){maxN = ceiling(3*median.prior.size)}
      x <- n:maxN
      a = optimize(f=fn1,interval=c(1,maxbeta),x,n,median.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      if(verbose){
       p=dfn(alpha,beta,x,n,effective.prior.df);
       cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
      }
      if(is.null(maxN.set)){
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
       abs(p[length(p)]/max(p,na.rm=TRUE) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p,na.rm=TRUE) > 0.01)+1]))
        x <- n:maxN
        a = optimize(f=fn1,interval=c(1,maxbeta),x,n,median.prior.size,
                     effective.prior.df,alpha,tol=0.01)
        beta <- a$minimum
        if(verbose){
       cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
        }
      }}
      if(is.null(maxN.set)){maxN = min(maxNmax,maxN)}
     }
     if(!is.null(mean.prior.size)){
      if(mean.prior.size < n){mean.prior.size = n}
      beta <- max(1.1,mean.prior.size/n - 1)
      if(is.null(alpha)) alpha=mean.prior.size/(mean.prior.size-n)
      beta <- 0.5*(1+alpha)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn2 <- function(beta,x,n,mean.prior.size,effective.prior.df,alpha){
       priorm <- dfn(alpha,beta,x,n,effective.prior.df)
       abs(mean.prior.size - sum(x*priorm)/sum(priorm,na.rm=TRUE))
      }
      if(is.null(maxN)){
      if(is.null(maxN.set)){maxN = ceiling(3*mean.prior.size)}
      x <- n:maxN
      a = optimize(f=fn2,interval=c(1,maxbeta),x,n,mean.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      if(verbose){
       p=dfn(alpha,beta,x,n,effective.prior.df);
       cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
      }
      if(is.null(maxN.set)){
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
       abs(p[length(p)]/max(p,na.rm=TRUE) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p,na.rm=TRUE) > 0.01)+1]))
        x <- n:maxN
        a = optimize(f=fn2,interval=c(1,maxbeta),x,n,mean.prior.size,
                    effective.prior.df,alpha,tol=0.01)
        beta <- a$minimum
        if(verbose){
       cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
        }
      }}
      if(is.null(maxN.set)){maxN = min(maxNmax,maxN)}
      }else{
      x <- 0:(maxN-1-n) + n
      a = optimize(f=fn2,interval=c(1,maxbeta),x,n,mean.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      }
     }
     if(!is.null(mode.prior.size)){
      if(mode.prior.size < n){mode.prior.size = n}
      beta <- 2*mode.prior.size/n - 1
      if(is.null(alpha)) alpha=mode.prior.size/(mode.prior.size-n)
      if(is.null(maxN)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
      fn3 <- function(beta,x,n,mode.prior.size,effective.prior.df,alpha){
       priorm <- dfn(alpha,beta,x,n,effective.prior.df)
       abs(mode.prior.size - (x+0.5)[which.max(priorm)])
      }
      if(is.null(maxN.set)){maxN = ceiling(3*mode.prior.size)}
      x <- n:maxN
      a = optimize(f=fn3,interval=c(1,maxbeta),x,n,mode.prior.size,
                   effective.prior.df,alpha,tol=0.01)
      beta <- a$minimum
      if(verbose){
       p=dfn(alpha,beta,x,n,effective.prior.df);
       cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
      }
      if(is.null(maxN.set)){
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
       abs(p[length(p)]/max(p,na.rm=TRUE) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p,na.rm=TRUE) > 0.01)+1]))
        x <- n:maxN
        a = optimize(f=fn3,interval=c(1,maxbeta),x,n,mode.prior.size,
                     effective.prior.df,alpha,tol=0.01)
        beta <- a$minimum
        if(verbose){
          cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
        }
      }}
      if(is.null(maxN.set)){maxN = min(maxNmax,maxN)}
     }
     if(!is.null(quartiles.prior.size)){
      if(quartiles.prior.size[1] < n){quartiles.prior.size[1] = n}
      if(quartiles.prior.size[2] < quartiles.prior.size[1]){
        aaa <- quartiles.prior.size[2]
        quartiles.prior.size[2] = quartiles.prior.size[1]
        quartiles.prior.size[1] = aaa
      }
      fn4 <- function(p,x,n,quartiles.prior.size,effective.prior.df){
       priorm <- dfn(exp(p[1]),exp(p[2]),x,n,effective.prior.df)
       sqrt((quartiles.prior.size[1] - (x+0.5)[match(TRUE,cumsum(priorm) >= 0.25)])^2+ 
            (quartiles.prior.size[2] - (x+0.5)[match(TRUE,cumsum(priorm) >= 0.75)])^2)
      }
      if(is.null(maxN.set)){maxN = ceiling(10*quartiles.prior.size[2])}
      x <- n:maxN
      a = optim(par=log(c(1,10)),fn=fn4,
        x=x,n=n,quartiles.prior.size=quartiles.prior.size,
        effective.prior.df=effective.prior.df,
        control=list(abstol=10))
      alpha <- exp(a$par[1])
      beta  <- exp(a$par[2])
      if(verbose){
       p=dfn(alpha,beta,x,n,effective.prior.df);
       cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
      }
#cat(sprintf("alpha=%f beta = %f value = %f maxN=%f E1=%f E3=%f Q1 = %f Q2 = %f\n",alpha,beta,a$value,maxN,(x+0.5)[match(TRUE,cumsum(p) >= 0.25)],(x+0.5)[match(TRUE,cumsum(p) >= 0.75)],quartiles.prior.size[1],quartiles.prior.size[2]))
      if(is.null(maxN.set)){
      while( {
       p=dfn(alpha,beta,x,n,effective.prior.df);
       abs(p[length(p)]/max(p,na.rm=TRUE) - 0.01)>0.005}){
        maxN <- round(maxN*(c(0.9,1.1)[(p[length(p)]/max(p,na.rm=TRUE) > 0.01)+1]))
        x <- n:maxN
        a = optim(par=a$par,fn=fn4,
         x=x,n=n,quartiles.prior.size=quartiles.prior.size,effective.prior.df=effective.prior.df,
         control=list(abstol=10))
        alpha <- exp(a$par[1])
        beta  <- exp(a$par[2])
        if(verbose){
       cat(sprintf("maxN= %d Est. Q1=%d Est. Q3=%d Q1 = %d Q2 = %d\n",maxN,round((x+0.5)[match(TRUE,cumsum(p) >= 0.25)]),round((x+0.5)[match(TRUE,cumsum(p) >= 0.75)]),quartiles.prior.size[1],quartiles.prior.size[2]))
        }
#cat(sprintf("alpha=%f beta = %f value = %f maxN=%f E1=%f E3=%f Q1 = %f Q2 = %f\n",alpha,beta,a$value,maxN,(x+0.5)[match(TRUE,cumsum(p) >= 0.25)],(x+0.5)[match(TRUE,cumsum(p) >= 0.75)],quartiles.prior.size[1],quartiles.prior.size[2]))
      }}
      if(is.null(maxN.set)){maxN = min(maxNmax,maxN)}
     }
     }
     cat(sprintf("Final alpha: %f, beta = %f max=%f\n",alpha, beta, maxbeta))
     if(is.null(alpha)) alpha=1
     if(is.null(beta)){
       warning("No prior information about the population size was specified! Using a prior mode of twice the sample size. Please specify prior information!", call. = FALSE)
       beta <- 3
     }
     if(is.infinite(beta) | is.na(beta)){
       stop("The sample size is incompatible with the specified prior information about the population size.", call. = FALSE)
     }
     if(is.null(maxN)&is.null(maxN.set)){maxN <- min(maxNmax,ceiling( n/(1-0.90^(1/beta)) ))}
     if(!is.null(maxN.set)){maxN <- maxN.set}
     if(is.null(N)){N <- min(maxNmax,ceiling( n/(1-0.5^(1/beta)) ))}
     x <- n:maxN
     lpriorm=dfn(alpha,beta,x,n,effective.prior.df);
     if(log){
      lpriorm <- log(lpriorm)
     }
     lpriorm
     },
    supplied={
     maxN <- supplied$maxN
     x <- n:maxN
     out <- supplied$sample
     outN <- out[,"N"]
     a=bgk_kde(outN,n=2^(ceiling(log(maxN-n)/log(2))),MIN=n,MAX=maxN)
     # Use an interpolating cubic spline
     posdensN <- spline(x=a[1,],y=a[2,],xout=x)$y
#    a=locfit( ~ lp(outN,nn=0.5))
#    posdensN <- predict(a, newdata=x)
     posdensN <- posdensN / sum(posdensN)
     lpriorm <- log(posdensN)
     lpriorm
     }
    )
#    End of switch to compute lpriorm
     x <- n:maxN
     if(log){
      priorm <- exp(lpriorm)
     }else{
      priorm <- lpriorm
     }
     if(priorsizedistribution!="flat"){
      mode.prior.size <- (x+0.5)[which.max(priorm)]
      mean.prior.size <- sum((x+0.5)*priorm)
      median.prior.size <- (x+0.5)[match(TRUE,cumsum(priorm) >= 0.5)]
      quartiles.prior.size[1] <- (x+0.5)[match(TRUE,cumsum(priorm) >= 0.25)]
      quartiles.prior.size[2] <- (x+0.5)[match(TRUE,cumsum(priorm) >= 0.75)]
     }else{
      mode.prior.size <- (maxN+n)/2
      mean.prior.size <- (maxN+n)/2
      median.prior.size <- (maxN+n)/2
      quartiles.prior.size[1] <-   (maxN+n)/4
      quartiles.prior.size[2] <- 3*(maxN+n)/4
     }
    if(is.null(N)){N <- mean(x)}
    if(verbose){cat(paste("The maximum prior population size is",maxN,"\n"))}
    out <- list(x=x,lprior=lpriorm,N=N,maxN=maxN,
         median.prior.size=median.prior.size,
         mean.prior.size=mean.prior.size,
         mode.prior.size=mode.prior.size,
         quartiles.prior.size=quartiles.prior.size,
	 mode.prior.sample.proportion=mode.prior.sample.proportion,
	 median.prior.sample.proportion=median.prior.sample.proportion,
	 alpha=alpha,beta=beta,
	 effective.prior.df=effective.prior.df,
         type=type)
    class(out) <- "sspse"
    out
}
