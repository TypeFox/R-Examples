#' @title Beta-Binomial distribution
#' @description Return the Beta-Binomial distribution with parameters m, csi and phi.
#' @aliases betar
#' @usage betar(m, csi, phi)
#' @param m Number of ordinal categories
#' @param csi  Feeling parameter of the Beta-Binomial distribution
#' @param phi Overdispersion parameter of the Beta-Binomial distribution 
#' @export betar
#' @return The vector of length \eqn{m} of the  Beta-Binomial distribution
#' @seealso   \code{\link{betabinomial}}
#' @references Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
#' \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786
#' @keywords distribution
#' @examples 
#' m<-9
#' csi<-0.8
#' phi<-0.2
#' pr<-betar(m, csi, phi)
#' plot(1:m,pr,type="h", main="Beta-Binomial distribution",xlab="Ordinal categories")
#' points(1:m,pr,pch=19)


betar <-function(m,csi,phi){
  betar<-rep(NA,m)
  km<-0:(m-2)
  betar[1]<-prod(1-(1-csi)/(1+phi*km))
  for(r in 1:(m-1)){
    betar[r+1]<-betar[r]*((m-r)/r)*((1-csi+phi*(r-1))/(csi+phi*(m-r-1)))
  }
  return(betar)
}

