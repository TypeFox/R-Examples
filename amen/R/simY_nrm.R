#' Simulate a normal relational matrix
#' 
#' Simulates a normal relational matrix
#' 
#' 
#' @usage simY_nrm(EY, rho, s2)
#' @param EY square matrix giving the expected value of the relational matrix
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @return a square matrix
#' @author Peter Hoff
#' @export simY_nrm
simY_nrm <-
function(EY,rho,s2) 
{
  YS<-simZ(EY,rho,s2) 
  diag(YS)<-NA 
  YS
}
