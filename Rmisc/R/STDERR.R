#' Standard Error
#' 
#' Calculates the standard error interval of a vector of data
#' 
#' @keywords univar
#' 
#' @param x a vector of data.
#' 
#' @return
#' \item{upper}{Upper bound of interval.}
#' \item{mean}{Mean of data.}
#' \item{lower}{Lower bound of interval.}
#' 
#' @export
#' 
#' @examples
#' STDERR(rnorm(100))
#' 
STDERR <-
function(x) {
  a<-mean(x)
  s<-sd(x)
  n<-length(x)
  error<-(s/sqrt(n))
  return(c(upper=a+error,mean=a,lower=a-error))
}
