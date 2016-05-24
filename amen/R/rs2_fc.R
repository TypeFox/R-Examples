#' Gibbs update for dyadic variance
#' 
#' Gibbs update for dyadic variance
#' 
#' 
#' @usage rs2_fc(E, rho)
#' @param E square residual relational matrix
#' @param rho current value of rho
#' @return a new value of s2
#' @author Peter Hoff
#' @export rs2_fc
rs2_fc <-
function(E,rho)
{
    H<-mhalf( solve(matrix(c(1,rho,rho,1),2,2)) )
    EM<-cbind(E[upper.tri(E)],t(E)[upper.tri(E)] ) %*%H
    1/rgamma(1, (length(EM)+1)/2 , (sum(EM^2)+1)/2 )
}
