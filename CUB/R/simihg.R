#' @title Simulation routine for IHG models
#' @aliases simihg
#' @description Generate n pseudo-random numbers according to the IHG distribution with given
#'  preference parameter \eqn{\theta}.
#' @keywords distribution
#' @usage simihg(n, m, theta)
#' @export simihg
#' @param n Number of simulated observations
#' @param m Number of ordinal categories
#' @param theta Preference parameter
#' @seealso  \code{\link{probihg}}
#' @examples
#' n<-300
#' m<-9
#' theta<-0.4
#' simulation<-simihg(n, m, theta)
#' plot(table(simulation), xlab="Number of categories", ylab="Frequencies")



simihg <-
function(n,m,theta){
  B<-(m-1)*theta/(1-theta)
  psi<-1-runif(n)^(1/B)
  vett<-1+rbinom(n,m-1,psi)
  return(vett)
}
