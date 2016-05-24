#' @title Probability distribution of an IHG model
#' @aliases probihg
#' @description Compute the probability distribution of an IHG model without covariates.
#' @keywords distribution
#' @export probihg
#' @usage probihg(m, theta)
#' @param m Number of ordinal categories
#' @param theta Preference parameter
#' @return The vector of the probability distribution of an IHG model
#' @seealso \code{\link{IHG}}
#' @references 
#' D'Elia A. (2003). Modelling ranks using the inverse hypergeometric distribution, 
#' \emph{Statistical Modelling: an International Journal}, \bold{3}, 65--78
#' @examples
#' m<-10
#' theta<-0.30
#' pr<-probihg(m, theta)
#' plot(1:m,pr,type="h",xlab="Ordinal categories")
#' points(1:m,pr,pch=19)


probihg <-
function(m,theta){
  pr<-rep(NA,m)
  pr[1]<-theta
  for(j in 1:(m-1)){
    pr[j+1]<-pr[j]*(1-theta)*(m-j)/(m-j-1+j*theta)
  }
  return(pr)
}
