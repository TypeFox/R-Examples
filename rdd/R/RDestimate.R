#' Regression Discontinuity Estimation
#' 
#' \code{RDestimate} supports both sharp and fuzzy RDD
#' utilizing the \pkg{AER} package for 2SLS regression
#' under the fuzzy design. Local linear regressions are performed
#' to either side of the cutpoint using the Imbens-Kalyanaraman
#' optimal bandwidth calculation, \code{\link{IKbandwidth}}.
#' 
#' @param formula the formula of the RDD. This is supplied in the
#' format of \code{y ~ x} for a simple sharp RDD, or \code{y ~ x | c1 + c2}
#' for a sharp RDD with two covariates. Fuzzy RDD may be specified as
#' \code{y ~ x + z} where \code{x} is the running variable, and 
#' \code{z} is the endogenous treatment variable. Covariates are then included in the 
#' same manner as in a sharp RDD.
#' @param data an optional data frame
#' @param subset an optional vector specifying a subset of observations to be used
#' @param cutpoint the cutpoint. If omitted, it is assumed to be 0.
#' @param bw a numeric vector specifying the bandwidths at which to estimate the RD. 
#' If omitted, the bandwidth is calculated using the Imbens-Kalyanaraman method, and then estimated
#' with that bandwidth, half that bandwidth, and twice that bandwidth. If only a single value is passed into the function,
#' the RD will similarly be estimated at that bandwidth, half that bandwidth, and twice that bandwidth.
#' @param kernel a string specifying the kernel to be used in the local linear fitting. 
#' \code{"triangular"} kernel is the default and is the "correct" theoretical kernel to be used for 
#' edge estimation as in RDD (Lee and Lemieux 2010). Other options are \code{"rectangular"}, 
#' \code{"epanechnikov"}, \code{"quartic"}, 
#' \code{"triweight"}, \code{"tricube"}, \code{"gaussian"} and \code{"cosine"}.
#' @param se.type this specifies the robust SE calculation method to use. Options are,
#' as in \code{\link{vcovHC}}, \code{"HC3"}, \code{"const"}, \code{"HC"}, \code{"HC0"}, 
#' \code{"HC1"}, \code{"HC2"}, \code{"HC4"}, \code{"HC4m"}, \code{"HC5"}. This option 
#' is overriden by \code{cluster}.
#' @param cluster an optional vector specifying clusters within which the errors are assumed
#' to be correlated. This will result in reporting cluster robust SEs. This option overrides
#' anything specified in \code{se.type}. It is suggested that data with a discrete running 
#' variable be clustered by each unique value of the running variable (Lee and Card 2008).
#' @param verbose will provide some additional information printed to the terminal.
#' @param frame logical. If \code{TRUE}, the data frame used in model fitting will be returned.
#' @param model logical. If \code{TRUE}, the model object will be returned.
#' @return \code{RDestimate} returns an object of \link{class} "\code{RD}".
#' The functions \code{summary} and \code{plot} are used to obtain and print a summary and plot of
#' the estimated regression discontinuity. The object of class \code{RD} is a list 
#' containing the following components:
#' \item{type}{a string denoting either \code{"sharp"} or \code{"fuzzy"} RDD.}
#' \item{est}{numeric vector of the estimate of the discontinuity in the outcome under a sharp design, 
#' or the Wald estimator in the fuzzy design for each corresponding bandwidth}
#' \item{se}{numeric vector of the standard error for each corresponding bandwidth}
#' \item{z}{numeric vector of the z statistic for each corresponding bandwidth}
#' \item{p}{numeric vector of the p value for each corresponding bandwidth}
#' \item{ci}{the matrix of the 95% confidence interval, \code{c("CI Lower Bound","CI Upper Bound")} 
#' for each corresponding bandwidth}
#' \item{bw}{numeric vector of each bandwidth used in estimation}
#' \item{obs}{vector of the number of observations within the corresponding bandwidth}
#' \item{call}{the matched call}
#' \item{na.action}{the observations removed from fitting due to missingness}
#' \item{model}{(if requested) For a sharp design, a list of the \code{lm} objects is returned.
#' For a fuzzy design, a list of lists is returned, each with two elements: \code{firststage}, the first stage \code{lm}
#' object, and \code{iv}, the \code{ivreg} object. A model is returned for each corresponding bandwidth.}
#' \item{frame}{(if requested) Returns the model frame used in fitting.}
#' @seealso \code{\link{summary.RD}}, \code{\link{plot.RD}}, \code{\link{DCdensity}} \code{\link{IKbandwidth}}, \code{\link{kernelwts}}, \code{\link{vcovHC}}, 
#' \code{\link{ivreg}}, \code{\link{lm}}
#' @details Covariates are problematic for inclusion in the regression
#' discontinuity design. This package allows their inclusion, but cautions
#' against them insomuch as is possible. When covariates are included in the
#' specification, they are simply included as exogenous regressors. In the
#' sharp design, this means they are simply added into the regression equation,
#' uninteracted with treatment. Likewise for the fuzzy design, in which they
#' are added as regressors in both stages of estimation.
#' @references Lee, David and Thomas Lemieux. (2010) "Regression Discontinuity Designs in Economics," \emph{Journal of Economic Literature}. 48(2): 281-355. \url{http://www.aeaweb.org/articles.php?doi=10.1257/jel.48.2.281}
#' @references Imbens, Guido and Thomas Lemieux. (2010) "Regression discontinuity designs: A guide to practice," \emph{Journal of Econometrics}. 142(2): 615-635. \url{http://dx.doi.org/10.1016/j.jeconom.2007.05.001}
#' @references Lee, David and David Card. (2010) "Regression discontinuity inference with specification error," \emph{Journal of Econometrics}. 142(2): 655-674. \url{http://dx.doi.org/10.1016/j.jeconom.2007.05.003}
#' @references Angrist, Joshua and Jorn-Steffen Pischke. (2009) \emph{Mostly Harmless Econometrics}. Princeton: Princeton University Press.
#' @import Formula AER sandwich lmtest
#' @importFrom stats model.frame na.pass complete.cases lm coef pnorm qnorm
#' as.formula
#' @include IKbandwidth.R 
#' @include kernelwts.R
#' @export
#' @author Drew Dimmery <\email{drewd@@nyu.edu}>
#' @examples
#' x<-runif(1000,-1,1)
#' cov<-rnorm(1000)
#' y<-3+2*x+3*cov+10*(x>=0)+rnorm(1000)
#' RDestimate(y~x)
#' # Efficiency gains can be made by including covariates
#' RDestimate(y~x|cov)

RDestimate<-function(formula, data, subset=NULL, cutpoint=NULL, bw=NULL, kernel="triangular", se.type="HC1", cluster=NULL, verbose=FALSE, model=FALSE, frame=FALSE) {
  call<-match.call()
  if(missing(data)) data<-environment(formula)
  formula<-as.Formula(formula)
  X<-model.frame(formula,rhs=1,lhs=0,data=data,na.action=na.pass)[[1]]
  Y<-model.frame(formula,rhs=0,lhs=NULL,data=data,na.action=na.pass)[[1]]
  if(!is.null(subset)){
    X<-X[subset]
    Y<-Y[subset]
    if(!is.null(cluster)) cluster<-cluster[subset]
  }
  if (!is.null(cluster)) {
    cluster<-as.character(cluster)
    robust.se <- function(model, cluster){
      M <- length(unique(cluster))
      N <- length(cluster)           
      K <- model$rank
      dfc <- (M/(M - 1)) * ((N - 1)/(N - K))
      uj  <- apply(estfun(model), 2, function(x) tapply(x, cluster, sum));
      rcse.cov <- dfc * sandwich(model, meat. = crossprod(uj)/N)
      rcse.se <- coeftest(model, rcse.cov)
      return(rcse.se[2,2])
    }
  }
  na.ok<-complete.cases(X)&complete.cases(Y)
  if(length(all.vars(formula(formula,rhs=1,lhs=F)))>1){
   type<-"fuzzy" 
   Z<-model.frame(formula,rhs=1,lhs=0,data=data,na.action=na.pass)[[2]]
   if(!is.null(subset)) Z<-Z[subset]
   na.ok<-na.ok&complete.cases(Z)
   if(length(all.vars(formula(formula,rhs=1,lhs=F)))>2)
     stop("Invalid formula. Read ?RDestimate for proper syntax")
  } else {
   type="sharp" 
  }
  covs<-NULL
  if(length(formula)[2]>1){
    covs<-model.frame(formula,rhs=2,lhs=0,data=data,na.action=na.pass)
    if(!is.null(subset)) covs<-subset(covs,subset)
    na.ok<-na.ok&complete.cases(covs)
    covs<-subset(covs,na.ok)
  }
  X<-X[na.ok]
  Y<-Y[na.ok]
  if(type=="fuzzy") Z<-as.double(Z[na.ok])
  
  if(is.null(cutpoint)) {
    cutpoint<-0
    if(verbose) cat("No cutpoint provided. Using default cutpoint of zero.\n")
  }
  if(frame) {
    if(type=="sharp") {
      if (!is.null(covs))
        dat.out<-data.frame(X,Y,covs)
      else
        dat.out<-data.frame(X,Y)
    } else {
      if (!is.null(covs))
        dat.out<-data.frame(X,Y,Z,covs)
      else
        dat.out<-data.frame(X,Y,Z)
    }
  }
  if(is.null(bw)) {
    bw<-IKbandwidth(X=X,Y=Y,cutpoint=cutpoint,kernel=kernel, verbose=verbose)
    bws<-c(bw,.5*bw,2*bw)
    names(bws)<-c("LATE","Half-BW","Double-BW")
  } else if (length(bw)==1) {
    bws<-c(bw,.5*bw,2*bw)
    names(bws)<-c("LATE","Half-BW","Double-BW")
  } else {
    bws<-bw
  }
  
  #Setup values to be returned
  o<-list()
  o$type<-type
  o$call<-call
  o$est<-vector(length=length(bws),mode="numeric")
  names(o$est)<-names(bws)
  o$bw<-as.vector(bws)
  o$se<-vector(mode="numeric")
  o$z<-vector(mode="numeric")
  o$p<-vector(mode="numeric")
  o$obs<-vector(mode="numeric")
  o$ci<-matrix(NA,nrow=length(bws),ncol=2)
  o$model<-list()
  if(type=="fuzzy") {
    o$model$firststage<-list()
    o$model$iv<-list()
  }
  o$frame<-list()
  o$na.action<-which(na.ok==FALSE)
  class(o)<-"RD"
  X<-X-cutpoint
  Xl<-(X<0)*X
  Xr<-(X>=0)*X
  Tr<-as.integer(X>=0)
  
  for(bw in bws){
    ibw<-which(bw==bws)
  #Subset to within the bandwidth, except for when using gaussian weighting
    sub<- X>=(-bw) & X<=(+bw)
  
    
  if(kernel=="gaussian") 
    sub<-TRUE

  w<-kernelwts(X,0,bw,kernel=kernel)
  o$obs[ibw]<-sum(w>0)
    
  if(type=="sharp"){
    if(verbose) {
      cat("Running Sharp RD\n")
      cat("Running variable:",all.vars(formula(formula,rhs=1,lhs=F))[1],"\n")
      cat("Outcome variable:",all.vars(formula(formula,rhs=F,lhs=1))[1],"\n")
      if(!is.null(covs)) cat("Covariates:",paste(names(covs),collapse=", "),"\n")
    }
    if(!is.null(covs)) {
      data<-data.frame(Y,Tr,Xl,Xr,covs,w)
      form<-as.formula(paste("Y~Tr+Xl+Xr+",paste(names(covs),collapse="+",sep=""),sep=""))
    } else {
      data<-data.frame(Y,Tr,Xl,Xr,w)
      form<-as.formula(Y~Tr+Xl+Xr)
    }

    mod<-lm(form,weights=w,data=subset(data,w>0))
    if(verbose==TRUE) {
      cat("Model:\n")
      print(summary(mod))
    }
    o$est[ibw]<-coef(mod)["Tr"]
    if (is.null(cluster)) {
      o$se[ibw]<-coeftest(mod,vcovHC(mod,type=se.type))[2,2]
    } else {
      o$se[ibw]<-robust.se(mod,cluster[na.ok][w>0])
    }
    o$z[ibw]<-o$est[ibw]/o$se[ibw]
    o$p[ibw]<-2*pnorm(abs(o$z[ibw]),lower.tail=F)
    o$ci[ibw,]<-c(o$est[ibw]-qnorm(.975)*o$se[ibw],o$est[ibw]+qnorm(.975)*o$se[ibw])

    if(model) o$model[[ibw]]=mod
    if(frame) o$frame[[ibw]]=dat.out

  } else {
    if(verbose){
      cat("Running Fuzzy RD\n")
      #CLEAN UP
      cat("Running variable:",all.vars(formula(formula,rhs=1,lhs=F))[1],"\n")
      cat("Outcome variable:",all.vars(formula(formula,rhs=F,lhs=1))[1],"\n")
      cat("Treatment variable:",all.vars(formula(formula,rhs=1,lhs=F))[2],"\n")
      if(!is.null(covs)) cat("Covariates:",paste(names(covs),collapse=", "),"\n")
    }

    if(!is.null(covs)) {
      data<-data.frame(Y,Tr,Xl,Xr,Z,covs,w)
      form<-as.Formula(paste(
              "Y~Z+Xl+Xr+",paste(names(covs),collapse="+"),
              "|Tr+Xl+Xr+",paste(names(covs),collapse="+"),sep=""))
      form1<-as.Formula(paste("Z~Tr+Xl+Xr+",paste(names(covs),collapse="+",sep="")))
    } else {
      data<-data.frame(Y,Tr,Xl,Xr,Z,w)
      form<-as.Formula(Y~Z+Xl+Xr|Tr+Xl+Xr)
      form1<-as.formula(Z~Tr+Xl+Xr)
    }
    
    mod1<-lm(form1,weights=w,data=subset(data,w>0))
    mod<-ivreg(form,weights=w,data=subset(data,w>0))
    if(verbose==TRUE) {
      cat("First stage:\n")
      print(summary(mod1))
      cat("IV-RD:\n")
      print(summary(mod))
    }
    o$est[ibw]<-coef(mod)["Z"]
    if (is.null(cluster)) {
      o$se[ibw]<-coeftest(mod,vcovHC(mod,type=se.type))[2,2]
    } else {
      o$se[ibw]<-robust.se(mod,cluster[na.ok][w>0])
    }
    o$z[ibw]<-o$est[ibw]/o$se[ibw]
    o$p[ibw]<-2*pnorm(abs(o$z[ibw]),lower.tail=F)
    o$ci[ibw,]<-c(o$est[ibw]-qnorm(.975)*o$se[ibw],o$est[ibw]+qnorm(.975)*o$se[ibw])
    
    if(model) {
      o$model$firststage[[ibw]]<-mod1
      o$model$iv[[ibw]]=mod
    }
    if(frame) o$frame=dat.out
  }
  }
  return(o)
}
