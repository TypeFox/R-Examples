#' Simulate missing values in a normal AME model
#' 
#' Simulates missing values of a sociomatrix under a normal AME model
#' 
#' 
#' @usage rZ_nrm_fc(Z, EZ, rho,s2, Y)
#' @param Z a square matrix, the current value of Z
#' @param EZ expected value of Z
#' @param rho dyadic correlation
#' @param s2 dyadic variance
#' @param Y square relational matrix
#' @return a square matrix, equal to Y at non-missing values
#' @author Peter Hoff
#' @export rZ_nrm_fc
rZ_nrm_fc<-function(Z,EZ,rho,s2,Y)
{
  ZS<-simY_nrm(EZ,rho,s2)
  diag(ZS)<-rnorm(nrow(Y),diag(EZ),sqrt(s2))
  Z[is.na(Y)]<-ZS[is.na(Y)]
  Z
}

