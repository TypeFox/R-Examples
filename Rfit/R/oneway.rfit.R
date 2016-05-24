#' Rank-based Oneway Analysis of Variance
#' 
#' Carries out a robust analysis of variance for a one factor design. Analysis
#' is based on the R estimates.
#' 
#' Carries out a robust one-way analysis of variance based on full model r fit.
#' 
#' @param y n by 1 response vector
#' @param g n by 1 vector representing group membership
#' @param scores an object of class 'scores'
#' @param p.adjust adjustment to the p-values, argument passed to p.adjust
#' @return \item{ fit }{ full model fit from rfit } \item{ est }{ Estimates }
#' \item{ se }{ Standard Errors } \item{ I }{ First Index } \item{ J }{ Second
#' Index } \item{p.value}{ p-values } \item{y}{ response vector } \item{g}{
#' vector denoting group membership }
#' @author Joseph McKean, John Kloke
#' @seealso \link{rfit}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @keywords anova nonparametric robust
#' @examples
#' 
#' 	data(quail)
#' 	oneway.rfit(quail$ldl,quail$treat)
#'  
#' @export oneway.rfit
oneway.rfit<-function(y,g,scores=Rfit::wscores,p.adjust='none') {

  ug<-unique(g)
  if( length(ug) < 3 ) stop('requires K > 2')
  K<-length(ug)
  nvec<-tapply(!is.na(y),g,sum)

  ### R fit ###
  x<-model.matrix(~as.factor(g)-1)
  x<-x[,2:ncol(x)]
  fit<-rfit(y~x)
  deltahat<-fit$betahat[-1]
  tauhat<-fit$tauhat

  ### R estimates of effect sizes ###
  kp<-length(deltahat)-1
  ind1<-rep(1:kp,times=seq(kp,1,by=-1))
  ind2<-ind1+sequence(seq(kp,1,by=-1))
  est<-c(deltahat,deltahat[ind1]-deltahat[ind2])

  ### Standard Errors ###
  kp<-length(deltahat)
  ind1<-rep(1:kp,times=seq(kp,1,by=-1))
  ind2<-ind1+sequence(seq(kp,1,by=-1))
  se<-tauhat*sqrt(1/nvec[ind1]+1/nvec[ind2])

  ### p-values ###
  pval<-p.adjust(2*pt(abs(est/se),length(y)-(kp+1),lower.tail=FALSE),
    method=p.adjust)
  # jk 20141119 - fix labels in p-value matrix 
  # pp<-matrix(nrow=K,ncol=K-1,dimnames=list(ug,ug[2:K]))
  # pp[lower.tri(pp)]<-pval
  pp<-matrix(nrow = K-1,ncol = K - 1, dimnames = list(ug[2:K],ug[1:(K-1)]))
  pp[lower.tri(pp,diag=TRUE)] <- pval
  DNAME <- paste(deparse(substitute(y)), "and", deparse(substitute(g)))
  ans<-list(method="Rfit",data.name=DNAME,p.value=pp,p.adjust.method=p.adjust)
  class(ans)<-"pairwise.htest"

  res<-list(fit=fit,est=est,se=se,I=ug[ind1],J=ug[ind2],p.value=pval,
    pp=ans,y=y,g=g)
  res$call<-match.call()
  class(res)<-"oneway.rfit"
  res

}
