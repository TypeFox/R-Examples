#' Simulate a network, i.e. a binary relational matrix
#' 
#' Simulates a network, i.e. a binary relational matrix
#' 
#' 
#' @usage simY_bin(EZ, rho)
#' @param EZ square matrix giving the expected value of the latent Z matrix
#' @param rho dyadic correlation
#' @return a square binary matrix
#' @author Peter Hoff
#' @export simY_bin
simY_bin <-
function(EZ,rho)
{
  ZS<-simZ(EZ,rho) 
  YS<-1*(ZS>0) ; diag(YS)<-NA
YS
}
