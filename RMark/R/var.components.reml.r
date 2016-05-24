#' Variance components estimation using REML or maximum likelihood
#' 
#' Computes estimated effects, standard errors and variance components for a
#' set of estimates
#' 
#' The function \code{\link{var.components}} uses method of moments to estimate
#' a single process variance but cannot fit a more complex example.  It can
#' only estimate an iid process variance.  However, if you have a more
#' complicated structure in which you have random year effects and want to
#' estimate a fixed age effect then \code{\link{var.components}} will not work
#' because it will assume an iid error rather than allowing a common error for
#' each year as well as an iid error.  This function uses restricted maximum
#' likelihood (reml) or maximum likelihood to fit a fixed effects model with an
#' optional random effects structure.  The example below provides an
#' illustration as to how this can be useful.
#' 
#' @param theta vector of parameter estimates
#' @param design design matrix for fixed effects combining parameter estimates
#' @param vcv estimated variance-covariance matrix for parameters
#' @param rdesign design matrix for random effect (do not use intercept form;
#' eg use ~-1+year instead of ~year); if NULL fits only iid error
#' @param initial initial values for variance components
#' @param interval interval bounds for log(sigma) to help optimization from
#' going awry
#' @param REML if TRUE uses reml else maximum likelihood
#' @return A list with the following elements \item{neglnl}{negative
#' log-likelihood for fitted model} \item{AICc}{small sample corrected AIC for
#' model selection} \item{sigma}{variance component estimates; if rdesign=NULL,
#' only an iid error; otherwise, iid error and random effect error}
#' \item{beta}{dataframe with estimates and standard errors of betas for
#' design} \item{vcv.beta}{variance-covariance matrix for beta}
#' @author Jeff Laake
#' @export
#' @examples
#' \donttest{
#' # This example is excluded from testing to reduce package check time
#' # Use dipper data with an age (0,1+)/time model for Phi
#' data(dipper)
#' dipper.proc=process.data(dipper,model="CJS")
#' dipper.ddl=make.design.data(dipper.proc,
#'    parameters=list(Phi=list(age.bins=c(0,.5,6))))
#' levels(dipper.ddl$Phi$age)=c("age0","age1+")
#' md=mark(dipper,model.parameters=list(Phi=list(formula=~time+age)))
#' # extract the estimates of Phi 
#' zz=get.real(md,"Phi",vcv=TRUE)
#' # assign age to use same intervals as these are not copied 
#' # across into the dataframe from get.real
#' zz$estimates$age=cut(zz$estimates$Age,c(0,.5,6),include=TRUE)
#' levels(zz$estimates$age)=c("age0","age1+")
#' z=zz$estimates
#' # Fit age fixed effects with random year component and an iid error
#' var.components.reml(z$estimate,design=model.matrix(~-1+age,z),
#'         zz$vcv,rdesign=model.matrix(~-1+time,z)) 
#' # Fitted model assuming no covariance structure to compare to 
#' # results with lme
#' xx=var.components.reml(z$estimate,design=model.matrix(~-1+age,z),
#'  matrix(0,nrow=nrow(zz$vcv),ncol=ncol(zz$vcv)),
#'  rdesign=model.matrix(~-1+time,z)) 
#' xx
#' sqrt(xx$sigmasq)
#' library(nlme)
#' lme(estimate~-1+age,data=z,random=~1|time)
#' }
var.components.reml=function(theta,design,vcv=NULL,rdesign=NULL,initial=NULL,interval=c(-25,10),REML=TRUE)
{
#  Arguments:
#     theta - vector of parameter estimates
#     design- design matrix for fixed effects
#     vcv   - estimated variance-covariance matrix for parameters
#     rdesign - design matrix for random effect (do not use intercept form; eg ~-1+year vs ~year)
#     initial - initial values for variance components
#     interval - interval bounds for log(sigma) to help optimization from going awry
#     REML - if TRUE uses reml else ml
#
#  Value: list with following elements
#     neglnl - negative log likelihood
#     AICc   - AICc value for model
#     sigma - process variance estimate
#     beta  - estimate of betas for design
#     vcv.beta - variance covariance matrix for beta
#
   if(is.null(initial))
     if(is.null(rdesign))
       initial=.1 
     else
       initial=c(.1,.1)
   if(nrow(design)!=length(theta)) stop("Number of rows of design matrix must match length of theta vector")
   if(!is.null(vcv))
   {
     if(nrow(vcv)!=length(theta)) stop("Number of rows of vcv matrix must match length of theta vector")
     if(ncol(vcv)!=length(theta)) stop("Number of columns of vcv matrix must match length of theta vector")
   }
   if(length(theta)<= ncol(design))stop("Length theta must exceed number of columns of design")
   if (!is.null(rdesign) && any(apply(rdesign, 2, FUN = function(x) all(abs(x - 1) < 1e-8))))
	   stop("Do not use intercept in formula for rdesign")
#
#  Reduce theta, design and vcv to the theta's that are used in the design
#   
   rn=(1:nrow(design))[apply(design,1,function(x)any(as.numeric(x)!=0))]
   theta=theta[rn]
   design=design[rn,,drop=FALSE]
   rdesign=rdesign[rn,,drop=FALSE]
   if(is.null(vcv))
     vcv=matrix(0,ncol=length(rn),nrow=length(rn))
   else
     vcv=vcv[rn,rn]
   n=length(rn)
   p=ncol(design)
# beta hat computation
   beta.hat = function(H, X, theta) return(solve(crossprod(X,solve(H,X)), crossprod(X,solve(H, theta))))
# Complete v-c matrix
   compute.H=function(alpha,vcv,rdesign)
   {
     sigmasq=exp(alpha)
     K=nrow(vcv)
     sigmamat=diag(sigmasq[1],nrow=K,ncol=K)
     if(!is.null(rdesign))
       sigmamat=sigmamat+rdesign%*%t(rdesign)*sigmasq[2]
     return(sigmamat+vcv)
   }
# reml negative log likelihood
   lnl.reml=function(alpha,vcv,X,theta,rdesign)
   {
      H=compute.H(alpha,vcv,rdesign)
      beta=beta.hat(H,X,theta)
      res=theta-X%*%beta
      lnl=.5*(determinant(H,logarithm=TRUE)$mod+crossprod(res,solve(H,res)))
      if(REML)lnl=lnl+.5*determinant(crossprod(X,solve(H,X)),logarithm=TRUE)$mod
      return(lnl)
   }
#  Find value of log(sigma) that minimizes reml negative log likelihood
   if(length(initial)==1)
   {
     soln=optimize(lnl.reml,interval=interval,vcv=vcv,theta=theta,X=design,rdesign=rdesign)
     estimate=soln$min
     neglnl=soln$objective
   }
   else
   {
     soln=optim(initial,lnl.reml,method="L-BFGS-B",lower=interval[1],upper=interval[2],vcv=vcv,theta=theta,X=design,rdesign=rdesign)
     estimate=soln$par
     neglnl=soln$value
   }
#  Compute beta-hat, std errors and v-c matrix at final values
   if(REML)
     neglnl=neglnl+.5*(n-p)*log(2*pi)
   else
     neglnl=neglnl+.5*n*log(2*pi)
   H=compute.H(estimate,vcv,rdesign)
   beta=beta.hat(H,design,theta)
   rownames(beta)=colnames(design)
   vcv.beta = solve(crossprod(design,solve(H,design)))
   rownames(vcv.beta)=colnames(design)
   colnames(vcv.beta)=colnames(design)
   beta=as.data.frame(beta)
   beta$se=sqrt(diag(vcv.beta))
   names(beta)=c("Estimate","SE")
   K=p+length(estimate)
   return(list(neglnl=as.vector(neglnl),AICc=as.vector(2*neglnl+2*K*(n/(n-K-1))),sigmasq=exp(estimate),beta=beta,vcv.beta=vcv.beta))
}
