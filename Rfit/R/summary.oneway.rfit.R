#' Provides a summary for the oneway anova based on an R fit.
#' 
#' Provides a summary for the oneway anova based on an R fit including a test
#' for main effects as tests for pairwise comparisons.
#' 
#' 
#' @param object an object of class 'oneway.rfit', usually, a result of a call
#' to 'oneway.rfit'
#' @param alpha Experimentwise Error Rate
#' @param method method used in confidence interval adjustment
#' @param \dots additional arguments
#' @author John Kloke, Joseph McKean
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' data(quail)
#' oneway.rfit(quail$ldl,quail$treat)
#' 
#' @export summary.oneway.rfit
summary.oneway.rfit<-function(object,alpha=0.05,method=confintadjust.methods,
  ...) {

  method<-match.arg(method)
  cv<-confintadjust(length(object$y),length(unique(object$g)),alpha,method)
  x<-cbind.data.frame(object$I,object$J,object$est,object$se,
    object$est-cv$cv*object$se,object$est+cv$cv*object$se)
  names(x)<-c('I','J','Estimate','St Err','Lower Bound CI', 'Upper Bound CI')
  res<-list(table=x,method=method)
  class(res)<-'summary.oneway.rfit'
  res

}
