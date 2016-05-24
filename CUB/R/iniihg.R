#' @title Moment estimate for the preference parameter of the IHG distribution
#' @description Compute the moment estimate of the preference parameter of the IHG distribution. 
#' This preliminary estimate is set as initial value within the optimization procedure for an IHG model
#'  fitting the observed frequencies.
#' @aliases iniihg
#' @usage iniihg(m, freq)
#' @param m Number of ordinal categories
#' @param freq Vector of length \eqn{m} of the absolute frequency distribution of the categories
#' @export iniihg
#' @return Moment estimator of the preference parameter \eqn{\theta}
#' @seealso \code{\link{IHG}}
#' @references D'Elia A. (2003). Modelling ranks using the inverse hypergeometric distribution, 
#' \emph{Statistical Modelling: An International Journal}, \bold{3}, 65--78
#' @keywords htest utilities
#' @examples
#' m<-9
#' freq<-c(70, 51, 48, 38, 29, 23, 12, 10, 5)
#' initheta<-iniihg(m, freq)



iniihg <-
function(m,freq){
  aver<-sum((1:m)*freq)/sum(freq)
  est<-(m-aver)/(1+(m-2)*aver) ### Moment estimator of theta
  return(est)
}
