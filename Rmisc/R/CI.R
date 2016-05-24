#' Confidence Interval
#' 
#' Calculates the confidence interval of a vector of data.
#' 
#' @keywords univar
#' 
#' @param x a vector of data
#' @param ci the confidence interval to be calculated
#' 
#' @return
#' \item{upper}{Upper bound of interval.}
#' \item{mean}{Mean of data.}
#' \item{lower}{Lower bound of interval.}
#' 
#' @export
#' 
#' @examples
#' CI(rnorm(100))
#' 
CI <-
function(x,ci=.95) {
  a<-mean(x)
  s<-sd(x)
  n<-length(x)
  error<-qt(ci+(1-ci)/2,df=n-1)*s/sqrt(n)
  return(c(upper=a+error,mean=a,lower=a-error))
}
