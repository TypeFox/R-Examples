#' @title Simulation routine for CUBE models
#' @aliases simcube
#' @description Generate \eqn{n} pseudo-random numbers according to the CUBE 
#' distribution with the given parameters.
#' @keywords distribution
#' @usage simcube(n, m, pai, csi, phi)
#' @export simcube
#' @param n Number of simulated observations
#' @param m Number of ordinal categories
#' @param pai Uncertainty parameter
#' @param csi Feeling parameter
#' @param phi Overdispersion parameter
#' @seealso \code{\link{probcube}}
#' @examples
#' n<-300
#' m<-9
#' pai<-0.7
#' csi<-0.4
#' phi<-0.1
#' simulation<-simcube(n,m,pai,csi,phi)
#' plot(table(simulation), xlab="Ordinal categories",ylab="Frequencies")


simcube <-
function(n,m,pai,csi,phi){
  prob<-probcube(m,pai,csi,phi)
  return(sample(1:m,n,prob,replace=T))
}
