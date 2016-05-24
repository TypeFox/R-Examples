# an example of match.arg 

confintadjust.methods<-c('none','bonferroni','tukey')



#' Confidence interval adjustment methods
#' 
#' Returns the critical value to be used in calculating adjusted confidence
#' intervals. Currently provides methods for Boneferroni and Tukey for
#' confidence interval adjustment methods as well as no adjustment.
#' 
#' Returns critial value based on one of the adjustment methods.
#' 
#' @aliases confintadjust confintadjust.methods
#' @param n sample size
#' @param k number of comparisons
#' @param alpha overall (experimentwise) type I error rate
#' @param method one of confintadjust.methods
#' @param \dots Additonal arguments.  Currently not used.
#' @return \item{cv}{critical value} \item{method}{the method used}
#' @author Joseph McKean, John Kloke
#' @references Hettmansperger, T.P. and McKean J.W. (2011), \emph{Robust
#' Nonparametric Statistical Methods, 2nd ed.}, New York: Chapman-Hall.
#' @examples
#' 
#' ##---- Should be DIRECTLY executable !! ----
#' ##-- ==>  Define data, use random,
#' ##--	or do  help(data=index)  for the standard data sets.
#' 
#' ## The function is currently defined as
#' function(n,k,alpha=0.05,method=confintadjust.methods,...) {
#' 	method<-match.arg(method)
#' 	cv<-switch(method, bonferroni = qt(1-alpha/choose(k,2),n-k),
#' 		tukey = qtukey(1-alpha,k,n-k)/sqrt(2),
#' 		none = qt(1-alpha/2,n-k)
#' 	)
#' 
#' 	res<-list(cv=cv,method=method)
#' 	res
#' 
#'   }
#' 
#' @export confintadjust
confintadjust<-function(n,k,alpha=0.05,method=confintadjust.methods,...) {
  method<-match.arg(method)
  cv<-switch(method, bonferroni = qt(1-alpha/choose(k,2),n-k),
    tukey = qtukey(1-alpha,k,n-k)/sqrt(2),
    none = qt(1-alpha/2,n-k)
  )

  res<-list(cv=cv,method=method)
  res

}


