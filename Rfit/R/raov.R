#raov<-function(f) {


#' R ANOVA
#' 
#' Returns full model fit and robust ANOVA table for all main effects and
#' interactions.
#' 
#' Based on reduction in dispersion tests for testing main effects and
#' interaction. Uses an algorithm described in Hocking (1985).
#' 
#' @param f an object of class formula
#' @param data an optional data frame
#' @param \dots additional arguments
#' @return \item{table}{Description of 'comp1'} \item{fit}{full model fit
#' returned from rfit} \item{residuals}{the residuals, i.e. y-yhat}
#' \item{fitted.values}{ yhat = x betahat} \item{call}{Call to the function}
#' @author Joseph McKean, John Kloke
#' @seealso \code{\link{rfit}}, \code{\link{oneway.rfit}}
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' 
#' Hocking, R. R. (1985), \emph{The Analysis of Linear Models}, Monterey,
#' California: Brooks/Cole.
#' @keywords anova nonparametric robust
#' @examples
#' 
#' raov(logSurv~Poison+Treatment,data=BoxCox)
#' 
#' @export raov
raov<-function(f,data=list(),...) {

# input formula : y ~ a + b +...+ z
# where a,b,...,z are k factors
	
  mf <- model.frame(formula = f,data=data)
#	X<-mf[,seq(ncol(mf),2,by=-1)]
#	X<-mf[,seq(2,ncol(mf))]
#	mf<-mf[do.call(order,X),]
	nf<-apply(mf[,2:ncol(mf)],2,function(foo) length(unique(foo)))

#  nf is number number of levels in each of the factors
#  mf is a data.frame of the form [y,a,b,...,z] sorted by z,...,b,a

  fit<-kwayr(nf,mf)

  # set up column names
  fnames<-names(mf)
  fnames<-fnames[2:length(fnames)]

  ind<-subsets(length(fnames))

  foo<-function(x,a) { paste(a[x==TRUE],collapse=':') }
  rownames(fit$tab)<-apply(ind,1,foo,a=fnames)

  colnames(fit$tab)<-c('DF','RD','Mean RD','F', 'p-value')
  res<-list(table=fit$tab,residuals=residuals(fit$fit),
    fitted.values=fitted.values(fit$fit),fit=fit$fit)

  res$call<-match.call()
  class(res)<-"raov"

  res

}
