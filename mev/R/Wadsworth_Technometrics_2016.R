#######################################################################################################
# Functions to create plots and tests in article
# Jennifer L. Wadsworth (2016), Technometrics
# "Exploiting structure of maximum likelihood estimators for extreme value threshold selection"
#
# Code by J.L. Wadsworth
#######################################################################################################

#' Wadsworth's univariate and bivariate exponential threshold diagnostics
#'
#' Function to produce diagnostic plots and test statistics for the
#' threshold diagnostics exploiting structure of maximum likelihood estimators
#' based on the non-homogeneous Poisson process likelihood
#'
#' @param xdat a numeric vector of data to be fitted.
#' @param model string specifying whether the univariate or bivariate diagnostic should be used. Either \code{nhpp}
#' for the univariate model, \code{exp} (\code{invexp}) for the bivariate exponential model with rate (inverse rate) parametrization. See details.
#' @param u optional; vector of candidate thresholds.
#' @param k number of thresholds to consider (if \code{u} unspecified).
#' @param q1 lowest quantile for the threshold sequence.
#' @param q2 upper quantile limit for the threshold sequence (\code{q2} itself is not used as a threshold,
#'  but rather the uppermost threshold will be at the \eqn{(q_2-1/k)}{q2-1/k} quantile).
#' @param par parameters of the NHPP likelihood. If \code{missing}, the \code{\link[evd]{fpot}} routine will be run to obtain values
#' @param M number of superpositions or "blocks" / "years" the process corresponds to (can affect the optimization)
#' @param nbs number of simulations used to assess the null distribution of the LRT, and produce the p-value
#' @param alpha significance level of the LRT
#' @param plots vector of strings indicating which plots to produce; \code{LRT}= likelihood ratio test, \code{WN} = white noise, \code{PS} = parameter stability
#' @param UseQuantiles logical; use quantiles as the thresholds in the plot?
#' @param pmar vector of length 4 giving the arguments for the plot margins in \code{par(mar=c(*,*,*,*))}.
#' @param tikz logical; if \code{TRUE}, axis labels are replaced with \code{LaTeX} code
#' @param ... additional parameters passed to \code{plot}.
#'
#' @details The function is a wrapper for the univariate (non-homogeneous Poisson process model) and bivariate exponential dependence model.
#' For the latter, the user can select either the rate or inverse rate parameter  (the inverse rate parametrization  works better for uniformity
#' of the p-value distribution under the \code{LR} test.
#'
#' There are two options for the bivariate diagnostic: either provide pairwise minimum of marginally
#' exponentially distributed margins or provide a \code{n} times 2 matrix with the original data, which
#' is transformed to exponential margins using the empirical distribution function.
#'
#' @references Wadsworth, J.L. (2016). Exploiting Structure of Maximum Likelihood Estimators for Extreme Value Threshold Selection, \emph{Technometrics}, \bold{58}(1), 116-126, \code{http://dx.doi.org/10.1080/00401706.2014.998345}.
#'
#' @author Jennifer L. Wadsworth
#' @return plots of the requested diagnostics and a list with components
#' \itemize{
#' \item \code{MLE}  maximum likelihood estimates from all thresholds
#' \item \code{Cov}  joint asymptotic covariance matrix for \eqn{\xi}{xi}, \eqn{\eta}{eta} or \eqn{\eta^{-1}}{1/eta}.
#' \item \code{WN}  values of the white noise process.
#' \item \code{LRT} values of the likelihood ratio test statistic vs threshold.
#' \item \code{pval} P-value of the likelihood ratio test.
#' \item \code{k}  final number of thresholds used.
#' \item \code{thresh} threshold selected by the likelihood ratio procedure.
#' \item \code{mle.u} maximum likelihood estimates from selected threshold.
#' }
#' @examples
#' \dontrun{
#' set.seed(123)
#' W.diag(rexp(1000), model="nhpp", k=30, q1=0)
#' # Parameter Stability only
#' W.diag(abs(rnorm(5000)), model="nhpp", k=30, q1=0, plots=c("PS"))
#' library(mvtnorm)
#' xbvn<-rmvnorm(6000, sigma=matrix(c(1,0.7,0.7,1),2,2))
#' # Transform margins to exponential manually
#' xbvn.exp<- -log(1-pnorm(xbvn))
#' W.diag(apply(xbvn.exp,1,min), model="exp", k=30, q1=0) #rate parametrization
#' W.diag(xbvn, model="exp", k=30, q1=0)
#' W.diag(apply(xbvn.exp,1,min), model="invexp", k=30, q1=0) #inverse rate parametrization
#' }
#' \dontrun{
#' library(ismev)
#' data(rain)
#' u <- quantile(rain, seq(0.85,0.99,by=0.01))
#' W.diag(xdat=rain, u=u, plots="PS")
#' }
#' @export
W.diag <- function(xdat, model=c("nhpp","exp","invexp"), u=NULL, k, q1=0, q2=1, par=NULL,
											M=NULL, nbs=1000, alpha=0.05, plots=c("LRT", "WN", "PS"),
											UseQuantiles=TRUE,pmar=c(5.5,7,3,3), tikz=FALSE, ...){
	model <- match.arg(model, c("nhpp","exp","invexp"))
	if(ncol(as.matrix(xdat)) == 2 && model %in% c("exp", "invexp")){
		xdat <- -log(1-apply(xdat, 2, function(y){rank(y, ties.method="average")/(length(y)+1)}))
		xdat <- pmin(xdat[,1],xdat[,2])
	}
	if(ncol(as.matrix(xdat))!=1){
		stop("Invalid input for `xdat`");
	}
	switch(model,
	nhpp=.NHPP.diag(xdat=xdat, u=u, k=k, q1=q1, q2=q2, par=par, M=M, nbs=nbs, alpha=alpha,
                      plots=plots, UseQuantiles=UseQuantiles, pmar=pmar,
                      tikz=tikz, ...),
	exp=.Expl.diag(x=xdat, u=u, k=k, q1=q1, q2=q2, nbs=nbs,
                      alpha=alpha, plots=plots,  UseQuantiles=UseQuantiles,
											param="Rate",pmar=pmar,...),
	invexp=.Expl.diag(x=xdat, u=u, k=k, q1=q1, q2=q2, nbs=nbs,
                      alpha=alpha, plots=plots,  UseQuantiles=UseQuantiles,
											param="InvRate",pmar=pmar,...)
	)
}


#' @importFrom ismev pp.fit
.NHPP.diag <- function(xdat, u=NULL, k, q1=0, q2=1, par=NULL, M=NULL, nbs=1000, alpha=0.05,
                      plots=c("LRT", "WN", "PS"), UseQuantiles=TRUE,pmar=c(5.5,7,3,3),
                      tikz=FALSE, ...){


	 	unull <- is.null(u)
		if(unull){
    	thresh <- quantile(xdat,q1)
     } else{
     	thresh <- min(u)
     }

  if(!unull){
  	k <- length(u)
  	}
	if(is.null(M)){
   M <- length(xdat[xdat>thresh])/3
   } #why M=nat/3 as default?
	 if(is.null(par)){
	{
  	  #require(ismev)
     	#Starting values for PP fit
#  	    uInd <- xdat > thresh
# 	    lrate <- sum(uInd)/length(xdat)
# 	    xdatu <- xdat[uInd]
# 	    in2 <- sqrt(6 * var(xdatu))/pi
# 	    in1 <- mean(xdatu) - 0.57722 * in2
# 	    in3 <- .Zhang_Stephens_posterior(xdatu-thresh)$mle['shape']
# 	    in2 <- exp(log(in2) + in3 * log(lrate))
# 	    in1 <- thresh - (in2/in3) * (lrate^(-in3) - 1)
# 			init.par <- list(loc=in1, scale=in2, shape=in3)
# 	    ppf <- evd::fpot(x=xdat, threshold=thresh, npp=length(xdat)/M, start=init.par, model="pp")
			ppf <- ismev::pp.fit(xdat=xdat,thresh=quantile(xdat,q1),npy=length(xdat)/M,show=F)
			par <- ppf$mle
#      par <- ppf$estimate
    }
	}
  par(mfrow=c(length(plots), 1), las=1, mar=pmar)

  J1 <- .Joint_MLE_NHPP(x=xdat, u=u, k=k, q1=q1, q2=q2, par=par,M=M)
  warn <- any(eigen(J1$Cov.xi)$val<=.Machine$double.eps)
  if(!unull && warn){
    stop("Estimated covariance matrix for xi not positive definite: try different thresholds")
  }

  while(any(eigen(J1$Cov.xi)$val<=.Machine$double.eps))
  {
    k <- k-1
    J1 <- .Joint_MLE_NHPP(x=xdat,k=k, q1=q1, q2=q2, par=par,M=M)
  }
  if(warn){
    warning(paste("Estimated covariance matrix for xi not positive definite for initial k. Final k:", k))
  }

  if(unull){
    u <- quantile(xdat,seq(q1,q2,len=k+1))
  }

  wn <- .C1(k)%*%J1$mle[,3]/sqrt(diag(.C1(k)%*%J1$Cov.xi%*%t(.C1(k))))
  nl <- .norm_LRT(x=wn,u=u[-c(1,k+1)])

  nlt <- NULL
  for(j in 1:nbs)
  {
    nlt[j] <- max(.norm_LRT(x=rnorm(k-1),u[-c(1,k+1)])[,2])
  }

  pval <- length(nlt[nlt>max(nl[,2])])/nbs

  if(pval<alpha){
    ustar <- nl[nl[,2]==max(nl[,2]),1]
  }  else{
    ustar <- min(u)
  }
  ind <- u[-(k+1)]==ustar
  theta.hat <- J1$mle[ind,]

  if(unull){
    qs <-  seq(q1, q2, len=k+1)[-(k+1)]
  }

  if(is.element("LRT",plots)){
    if(UseQuantiles  &&  unull){
      plot(qs, c(rep(NA,2),nl[,2]),xlab="Quantile",ylab="LR statistic", main=paste("p-value:", pval),...)
    } else{
      plot(u[-c(k+1)], c(rep(NA,2),nl[,2]),xlab="Threshold",ylab="LR statistic", main=paste("p-value:", pval),...)
    }
  }

  if(is.element("WN",plots)){
    if(UseQuantiles  &&  unull){
      plot(qs, c(NA,wn),xlab="Quantile",ylab="White noise",...)
      abline(h=0,col=2)
      abline(v=mean(xdat<=ustar),col=4)
    } else{
      plot(u[-c(k+1)], c(NA,wn),xlab="Threshold",ylab="White noise",...)
      abline(h=0,col=2)
      abline(v=ustar,col=4)
    }
  }

  if(is.element("PS",plots)){
    TradCI <- cbind(J1$mle[,3]-qnorm(0.975)*sqrt(diag(J1$Cov.xi)),J1$mle[,3]+qnorm(0.975)*sqrt(diag(J1$Cov.xi)))
    if(UseQuantiles  &&  unull){
      plot(qs,J1$mle[,3],ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(xi)),...)
      lines(qs,TradCI[,1], lty=2)
      lines(qs,TradCI[,2], lty=2)
      abline(v=mean(xdat<=ustar),col=4)
    }  else{
      plot(u[-(k+1)],J1$mle[,3],ylim=c(min(TradCI[,1]),max(TradCI[,2])),
           xlab="Threshold", ylab=ifelse(tikz,"$\\hat{\\xi}$",
                                         expression(hat(xi))),...)
      lines(u[-(k+1)],TradCI[,1], lty=2)
      lines(u[-(k+1)],TradCI[,2], lty=2)
      abline(v=ustar,col=4)}
  }

  return(list(MLE=J1$mle, Cov=J1$Cov.xi, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
}


#############################################################################################################

.Expl.diag <- function(x, u= NULL, k, q1, q2=1, nbs=1000,
                      alpha=0.05, plots=c("LRT", "WN", "PS"),
                      UseQuantiles=TRUE, param="InvRate",pmar=c(5.5,7,3,3),...){
  unull <- is.null(u)
  if(!unull){
    k <- length(u)
  }

  par(mfrow=c(length(plots), 1),mar=pmar, las=1)
  J1 <- .Joint_MLE_Expl(x=x,u=u,k=k, q1=q1, q2=q2, param=param)
  warn <- any(eigen(J1$Cov)$val<=.Machine$double.eps)
  if(!unull && warn){
    stop("Estimated covariance matrix for eta not positive definite: try different thresholds")
  }

  while(any(eigen(J1$Cov)$val<=.Machine$double.eps))  {
    k <- k-1
    J1 <- .Joint_MLE_Expl(x=x,k=k, q1=q1, q2=q2, param=param)
  }
  if(warn){
    warning(paste("Estimated covariance matrix for 1/eta not positive definite for initial k. Final k:", k))
  }
  if(unull){
    u <- quantile(x,seq(q1,1,len=k+1))
  }
  wn <- .C1(k)%*%J1$mle/sqrt(diag(.C1(k)%*%J1$Cov%*%t(.C1(k))))
  nl <- .norm_LRT(x=wn,u=u[-c(1,k+1)])

  nlt <- NULL
  for(j in 1:nbs){
    nlt[j] <- max(.norm_LRT(x=rnorm(k-1),u[-c(1,k+1)])[,2])
  }

  pval <- length(nlt[nlt>max(nl[,2])])/nbs

  if(pval<alpha){
    ustar <- nl[nl[,2]==max(nl[,2]),1]
  }  else{
    ustar <- min(u)
  }
  ind <- u[-(k+1)]==ustar
  theta.hat <- J1$mle[ind]

  if(unull){
    qs <-  seq(q1, q2, len=k+1)[-(k+1)]
  }

  if(is.element("LRT",plots)){
    if(unull && UseQuantiles){
      plot(qs, c(NA,NA,nl[,2]),xlab="Quantile",ylab="LR statistic", main=paste("p-value:", pval),...)
    } else{
      plot(u[-c(k+1)], c(NA,NA,nl[,2]),xlab="Threshold",ylab="LR statistic", main=paste("p-value:", pval),...)
    }
  }

  if(is.element("WN",plots)){
    if(unull && UseQuantiles){plot(qs, c(NA,wn),xlab="Quantile",ylab="White noise",...)
      abline(h=0,col=2)
      abline(v=mean(x<=ustar),col=4)
    } else{
      plot(u[-c(k+1)], c(NA,wn),xlab="Threshold",ylab="White noise",...)
      abline(h=0,col=2)
      abline(v=ustar,col=4)
    }
  }

  if(is.element("PS",plots)){
    TradCI <- cbind(J1$mle-qnorm(0.975)*sqrt(diag(J1$Cov)),J1$mle+qnorm(0.975)*sqrt(diag(J1$Cov)))
    if(UseQuantiles)    {
      if(param=="InvRate"){
        plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(eta)),...)
      } else if(param=="Rate"){
        plot(qs,J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Quantile", ylab=expression(hat(theta)),...)}
      lines(qs,TradCI[,1], lty=2)
      lines(qs,TradCI[,2], lty=2)
      abline(v=mean(x<=ustar),col=4)
    } else{
      if(param=="InvRate"){
        plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])),
             xlab="Threshold", ylab=expression(hat(eta)),...)
      } else if(param=="InvRate"){
        plot(u[-(k+1)],J1$mle,ylim=c(min(TradCI[,1]),max(TradCI[,2])), xlab="Threshold",
             ylab=expression(hat(theta)),...)
      }
      lines(u[-(k+1)],TradCI[,1], lty=2)
      lines(u[-(k+1)],TradCI[,2], lty=2)
      abline(v=ustar,col=4)
    }
  }

  return(list(MLE=J1$mle, Cov=J1$Cov, WN=wn, LRT = nl, pval = pval, k=k, thresh=ustar, mle.u=theta.hat))
}







#######################################################################################################

#' Joint maximum likelihood estimation for exponential model
#'
#'
#' Calculates the MLEs of the rate parameter, and joint asymptotic covariance matrix of these MLEs
#' over a range of thresholds as supplied by the user.
#'
#'@param x vector of data
#'@param u vector of thresholds. If not supplied, then \code{k}
#' thresholds between quantiles (\code{q1}, \code{q2}) will be used
#'@param k number of thresholds to consider if u not supplied
#'@param q1 lower quantile to consider for threshold
#'@param q2 upper quantile to consider for threshold
#'@param param character specifying \code{"InvRate"} or \code{"Rate"}
#' for either inverse rate parameter / rate parameter, respectively
#'
#' @author Jennifer L. Wadsworth
#'
# Value: (returns:)
#'@return a list with
#' \itemize{
#' \item mle vector of MLEs above the supplied thresholds
#' \item cov joint asymptotic covariance matrix of these MLEs
#' }
.Joint_MLE_Expl <- function(x, u=NULL, k, q1, q2=1, param){
  if(!is.element(param, c("InvRate","Rate"))){
    stop("param should be one of InvRate or Rate")
    }
  if(!is.null(u)){
    k <- length(u)
    x <- x[x>u[1]]
    # add threshold above all data, to "close" final region
    u <- c(u,max(x)+1)
  } else{
    u <- quantile(x,seq(q1,q2,len=k+1))
    }

  I <- n <- m <- thetahat <- NULL

  for(i in 1:k){
    if(param=="InvRate"){
      thetahat[i] <- mean(x[x>=u[i]]- u[i])
    } else if(param=="Rate"){
      thetahat[i] <- 1/mean(x[x>=u[i]]- u[i])
    }
    n[i] <- length(x[x>=u[i]&x<=u[i+1]])
  }
  for(i in 1:k) {
    m[i] <- sum(n[i:k])
    I[i] <-  1/thetahat[i]^2
  }

  Tcov <- matrix(0,k,k)
  for(i in 1:k) {
    for(j in 1:k) {
      Tcov[i,j] <- 1/(I[min(i,j)]*m[min(i,j)])
    }
  }
  CovT <- Tcov
  return(list(mle=thetahat, Cov=CovT))
}

#####################################################################################
#' Negative log-likelihood of the non-homogeneous Poisson process (NHPP)
#'
#' Negative log likelihood for the NHPP model, to minimize for MLEs
#'
#' @param theta parameter vector (\eqn{\mu}, \eqn{sigma}, \eqn{xi})
#' @param x data vector
#' @param u threshold
#' @param M number of superpositions or "blocks" / "years" the process corresponds to (affects estimation of \eqn{mu}, #' \eqn{sigma}, but these can be changed post-hoc to correspond to any number)
#'
#' @author Jennifer L. Wadsworth
#'
#' @return the value of the negative log-likelihood value
.nhpp_nll <- function(theta,x,u,M){
  x <- x[x>u]
  mu <- theta[1]; sig <- theta[2]; xi <- theta[3]
  Nu <- length(x)

  if(sig<=0|| any(1+xi*(x-mu)/sig < 0)){
    return(10e10)
  }  else{
    if(abs(xi)>1e-10)    {
      nll <-  (1/xi+1)*sum(log(1+xi*(x-mu)/sig)) + Nu*log(sig) + M*(1+xi*(u-mu)/sig)^(-1/xi)
    } else{
      nll <- Nu*log(sig)+sum((x-mu)/sig) + M*exp(-(u-mu)/sig)
    }
    return(nll)
  }
}



#####################################################################################


#' Joint maximum likelihood for the non-homogeneous Poisson Process
#'
#' Calculates the MLEs of the parameters (\eqn{\mu}, \eqn{\sigma}, \eqn{\xi}), and joint
#' asymptotic covariance matrix of these MLEs over a range of thresholds as supplied by the user.
#'@param x vector of data
#'@param u optional vector of thresholds. If not supplied, then k thresholds between quantiles (q1, q2) will be used
#'@param k number of thresholds to consider if \code{u} not supplied
#'@param q1 lower  quantile to consider for threshold
#'@param q2 upper quantile to consider for threshold. Default to 1
#'@param par starting values for the optimization
#'@param  M  number of superpositions or "blocks" / "years" the process corresponds to.
#' It affects the estimation of \eqn{mu} and \eqn{sigma},
#' but these can be changed post-hoc to correspond to any number)
#'
#' #' @author Jennifer L. Wadsworth
#'@return a list with components
#'\itemize{
#'\item mle matrix of MLEs above the supplied thresholds; columns are (\eqn{\mu}, \eqn{\sigma}, \eqn{\xi})
#'\item Cov.all joint asymptotic covariance matrix of all MLEs
#'\item Cov.mu joint asymptotic covariance matrix of MLEs for \eqn{\mu}
#'\item Cov.sig joint asymptotic covariance matrix of MLEs for \eqn{\sigma}
#'\item Cov.xi joint asymptotic covariance matrix of MLEs for \eqn{\xi}
#'}
.Joint_MLE_NHPP <- function(x, u=NULL, k, q1, q2=1, par, M){
  if(!is.null(u)) {
    k <- length(u)
    x <- x[x>u[1]]
    # add threshold above all data, to "close" final region
    u <- c(u,max(x)+1)
  }  else{
    u <- quantile(x,seq(q1,q2,len=k+1))
  }

  I <- Iinv <- list()
  thetahat <- matrix(NA,ncol=3,nrow=k)

  for(i in 1:k){
    opt <- optim(.nhpp_nll, par=par, x=x, u=u[i],M=M, hessian=F)
    thetahat[i,] <- opt$par

    ###
    # Deal with xi <- 0.5
    if(thetahat[i,3]>-0.5){
      ###
      I[[i]] <-  .E_Info_Mat(theta=opt$par, u=u[i], M=M)$EIM
      Iinv[[i]] <- solve(I[[i]])
    } else{
      I[[i]] <- Iinv[[i]] <- matrix(0,3,3)
    }
  }

  Wcov <- list()
  Wcov1 <- NULL
  for(i in 1:k){
    Wcov[[i]] <- matrix(0,3,3)
    for(j in 1:k)    {
      Wcov[[i]] <- cbind(Wcov[[i]],Iinv[[min(i,j)]])
    }
    Wcov1 <- rbind(Wcov1,Wcov[[i]])
  }
  Wcov1 <- Wcov1[,-c(1:3)]

  CovT <- Wcov1

  Cov.mu <-  CovT[seq(1,3*k,by=3),seq(1,3*k,by=3)]
  Cov.sig <- CovT[seq(2,3*k,by=3),seq(2,3*k,by=3)]
  Cov.xi <-  CovT[seq(3,3*k,by=3),seq(3,3*k,by=3)]

  return(list(mle=thetahat, Cov.all=CovT, Cov.mu=Cov.mu, Cov.sig = Cov.sig, Cov.xi=Cov.xi))
}


###################################################################################

# norm_LRT

# Details:

# Evaluates the likelihood ratio statistics for testing white noise

# Arguments:

# x - vector of white noise process (WNP, usually normalized estimates of \eqn{xi} or the exponential rate parameter \eqn{1/\eta})
# u - vector of thresholds that are associated to the WNP


.norm_LRT <- function(x,u){
  l <- length(u)
  v <- u[-c(1)] # means two or more obs available for std dev calculation
  lr <- NULL
  for(i in 1:length(v)){
    n1  <- length(x[u<=v[i]])
    num <- .nll_norm(theta=c(mean(x[u<=v[i]]),  sd(x[u<=v[i]])*sqrt((n1-1)/n1)), x=x[u<=v[i]])
    den <- .nll_norm(theta=c(0,1), x=x[u<=v[i]])
    lr[i] <- -2*(num-den)
  }
  return(cbind(v,lr))
}



###################################################################################

# nll_norm - negative log likelihood for the normal distribution

.nll_norm <- function(theta,x){
  if(theta[2]<0){
    return(10e10)
  }  else{
    return(-sum(dnorm(x,mean=theta[1],sd=theta[2],log=T)))
  }
}



###################################################################################

#' Contrast matrix
#'
#' Produces a contrast matrix with (1,-1) elements running down the two diagonals
#'
#'@param k number of columns (the number of rows is \code{k-1})
#'
#'@return a \code{k-1} x \code{k} contrast matrix
#'
.C1 <- function(k){
  C <- diag(x=1, nrow=k-1, ncol=k)
  C[row(C)+1 == col(C)] <- -1
  return(C)
}


###################################################################################

#' Functions for the expected information matrix used in .Joint_MLE_NHPP
#'
#' Functions with names of form \code{"d2ldmu2"} are derivatives (with expectation over the random NUMBER of points
#' already incorporated). Functions with names of form \code{"i_d2ldmu2"}, are in a form ready for integration;
#' this integration yields the expectation over the random LOCATIONS of the points.
#'
#' @param x dummy variable over which to integrate
#' @param mu location parameter
#' @param sig scale parameter
#' @param shape parameter
#' @param u threshold above which NHPP model assumed
#' @param M number of superpositions or "blocks" / "years" the process corresponds to
#'
.d2ldmu2 <- function(x, mu, sig, xi, u, M){
  c1 <- (1+(xi/sig)*(u-mu))^(-1/xi)
  p1 <- M*c1*(xi*(1+xi)/sig^2)*(1+(xi/sig)*(x-mu))^(-2)
  p2 <-  -M*((1+xi)/sig^2)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  return(p1+p2)
}

.i_d2ldmu2 <- function(x, mu, sig, xi, u, M){
  .d2ldmu2(x=x,mu=mu,sig=sig,xi=xi,u=u,M=M)* (1+(xi/sig)*(u-mu))^(1/xi) *
    (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



.d2ldmudsig <- function(x, mu, sig, xi, u, M){
  c1 <- (1+(xi/sig)*(u-mu))^(-1/xi)
  p1 <- M*c1*(xi*(1+xi)/sig^3)*(x-mu)*(1+(xi/sig)*(x-mu))^(-2)
  p2 <- M*c1*(-(1+xi)/sig^2)*(1+(xi/sig)*(x-mu))^(-1)
  p3 <- M*(1/sig^2)*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  p4 <- -M*((1+xi)/sig^3)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  return(p1+p2+p3+p4)
}

.i_d2ldmudsig <- function(x, mu, sig, xi, u, M){
  .d2ldmudsig(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) *
    (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}

#########

.d2ldmudxi <- function(x, mu, sig, xi, u, M){
  c1 <- (1+(xi/sig)*(u-mu))^(-1/xi)
  p1 <- M*c1* (1/sig)*(1+(xi/sig)*(x-mu))^(-1)
  p2 <- M*c1*(-(1+xi)/sig^2)*(x-mu)*(1+(xi/sig)*(x-mu))^(-2)
  p3.1 <- (1/sig)*((-1/xi^2)*log(1+(xi/sig)*(u-mu)) + (1/sig)*(1/xi+1)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))
  p3 <- M*p3.1*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  return(p1+p2+p3)
}

.i_d2ldmudxi <- function(x, mu, sig, xi, u, M){
  .d2ldmudxi(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) *
    (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



##########


.d2ldsig2 <- function(x, mu, sig, xi, u, M){
  c1 <- (1+(xi/sig)*(u-mu))^(-1/xi)
  p1 <- M*c1*(-2*(1+xi)/sig^3)*(x-mu)*(1+(xi/sig)*(x-mu))^(-1)
  p2 <- M*c1*(xi*(1+xi)/sig^4)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)
  p3 <-  M*(2/sig^3)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  p4 <-  -M*((1+xi)/sig^4)*((u-mu)^2)*(1+(xi/sig)*(u-mu))^(-1/xi-2)
  p5 <- M*c1*1/sig^2
  return(p1+p2+p3+p4+p5)
}

.i_d2ldsig2 <- function(x, mu, sig, xi, u, M){
  .d2ldsig2(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) *
    (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}



#########

.d2ldsigdxi <- function(x, mu, sig, xi, u, M){
  c1 <- (1+(xi/sig)*(u-mu))^(-1/xi)
  p1 <-  M*c1*((x-mu)/sig^2)*(1+(xi/sig)*(x-mu))^(-1)
  p2 <- M*c1*(-(1+xi)/sig^3)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)
  p3.1 <- ((u-mu)/sig^2)*((-1/xi^2)*log(1+(xi/sig)*(u-mu)) +
                            (1/sig)*(1/xi+1)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))
  p3 <- M*p3.1*(1+(xi/sig)*(u-mu))^(-1/xi-1)
  return(p1+p2+p3)
}


.i_d2ldsigdxi <- function(x, mu, sig, xi, u, M){
  .d2ldsigdxi(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)*
    (1+(xi/sig)*(u-mu))^(1/xi) * (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}


#########

.d2ldxi2 <- function(x, mu, sig, xi, u, M){
  c1 <- (1+(xi/sig)*(u-mu))^(-1/xi)
  p1 <-  M*c1*(-2/xi^3)*log(1+(xi/sig)*(x-mu))#
  p2 <-  M*c1*2*((x-mu)/(sig*xi^2))*(1+(xi/sig)*(x-mu))^(-1)#
  p3 <- M*c1*((1/xi+1)/sig^2)*((x-mu)^2)*(1+(xi/sig)*(x-mu))^(-2)#
  p4.1 <- ((-1/xi^2)*log(1+(xi/sig)* (u-mu)) + (1/sig)*(1/xi)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1))#
  p4 <- -(p4.1^2)*M*(1+(xi/sig)*(u-mu))^(-1/xi)#
  p5.1 <- ((2/xi^3)*log(1+(xi/sig)*(u-mu)) - (2/sig)*(1/xi^2)*(u-mu)*(1+(xi/sig)*(u-mu))^(-1)-
             (1/sig^2)*(1/xi)*((u-mu)^2)*(1+(xi/sig)*(u-mu))^(-2))#
  p5 <- (p5.1)*M*(1+(xi/sig)*(u-mu))^(-1/xi)#

  return(p1+p2+p3+p4+p5)
}

.i_d2ldxi2 <- function(x, mu, sig, xi, u, M){
  .d2ldxi2(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1+(xi/sig)*(u-mu))^(1/xi) *
    (1/sig)*(1+(xi/sig)*(x-mu))^(-1/xi-1)
}

## version of d2ldxi2 and i_d2ldxi2 with limit as xi \to 0. (This derivative has most trouble with small absolute
## values of xi.) Also possible to take such limits in other derivatives, but not implemented as they are generally
## less problematic.


.d2ldxi2_xi0 <- function(x, mu, sig, xi, u, M){
  c1 <- exp(-(u-mu)/sig)
  q1 <- ((x-mu)/sig)
  q2 <- ((u-mu)/sig)
  p1 <- -((2/3)*q1^3-q1^2)*M*c1
  p2 <- -(-(2/3)*q2^3+(1/4)*q2^4)*M*c1

  return(p1+p2)
}

.i_d2ldxi2_xi0 <- function(x, mu, sig, xi, u, M){
  .d2ldxi2_xi0(x=x,mu=mu,sig=sig,xi=xi,u=u, M=M)* (1/sig)*exp(-(x-u)/sig)
}



#############################################################################################
#' Numerically integrated expected information matrix of the NHPP
#'
#' Calculates the numerically-integrated expected information matrix for an NHPP with specified parameters
#' @param theta vector of parameters (\eqn{\mu}, \eqn{\sigma}, \eqn{\xi})
#' @param u threshold for NHPP
#' @param M number of superpositions or "blocks" / "years" the process corresponds to
#' @author Jennifer L. Wadsworth
#'
#' @return a list with components
#' \itemize{
#' \item EIM  expected information matrix
#' \item Errors  vector of errors from the numerical integration of the 6 unique components
#' }
.E_Info_Mat <- function(theta, u, M){
  if(theta[3]<0){
    up <- theta[1] - theta[2]/theta[3]
  } else {
    up <- Inf
  }
  Ed2ldmu2 <- integrate(.i_d2ldmu2,lower=u,upper=up,mu=theta[1], sig=theta[2],
                        xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldmudsig <- integrate(.i_d2ldmudsig,lower=u,upper=up,mu=theta[1], sig=theta[2],
                           xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldmudxi <- integrate(.i_d2ldmudxi,lower=u,upper=up,mu=theta[1], sig=theta[2],
                          xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldsig2 <- integrate(.i_d2ldsig2,lower=u,upper=up,mu=theta[1], sig=theta[2],
                         xi=theta[3], u=u,M=M, abs.tol=0)
  Ed2ldsigdxi <- integrate(.i_d2ldsigdxi,lower=u,upper=up,mu=theta[1], sig=theta[2],
                           xi=theta[3], u=u,M=M, abs.tol=0)
  if(abs(theta[3])>0.0001){
    Ed2ldxi2 <- integrate(.i_d2ldxi2,lower=u,upper=up,mu=theta[1], sig=theta[2],
                          xi=theta[3], u=u,M=M, abs.tol=0)
  } else{
    Ed2ldxi2 <- integrate(.i_d2ldxi2_xi0,lower=u,upper=up,mu=theta[1], sig=theta[2],
                          xi=theta[3], u=u,M=M, abs.tol=0)
  }

  Errors <- c(Ed2ldmu2$abs.err,Ed2ldmudsig$abs.err,Ed2ldmudxi$abs.err,Ed2ldmudsig$abs.err,Ed2ldsig2$abs.err,
              Ed2ldsigdxi$abs.err,Ed2ldmudxi$abs.err,Ed2ldsigdxi$abs.err,Ed2ldxi2$abs.err)

  EIM <- -matrix(c(Ed2ldmu2$val,Ed2ldmudsig$val,Ed2ldmudxi$val,Ed2ldmudsig$val,Ed2ldsig2$val,Ed2ldsigdxi$val,
                   Ed2ldmudxi$val,Ed2ldsigdxi$val,Ed2ldxi2$val), byrow=T, nrow=3)

  return(list(EIM=EIM,Errors=Errors))
}
