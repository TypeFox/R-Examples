#' Gibbs update for dyadic variance with independent replicate relational data
#' 
#' Gibbs update for dyadic variance with independent replicate relational data
#' 
#' 
#' @usage rs2_rep_fc(E.T, rho)
#' @param E.T Array of square residual relational matrix series. The third
#' dimension of the array is for different replicates. Each slice of the array
#' according to the third dimension is a square residual relational matrix
#' @param rho current value of rho
#' @return a new value of s2
#' @author Peter Hoff, Yanjun He
#' @export rs2_rep_fc
rs2_rep_fc <-
  function(E.T,rho)
  {
    N<-dim(E.T)[3]
    H<-mhalf( solve(matrix(c(1,rho,rho,1),2,2)) )
    EM<-NULL
    for (t in 1:N){
      E<-E.T[,,t]
      EM<-rbind(EM,cbind(E[upper.tri(E)],t(E)[upper.tri(E)] ) %*%H)
    } 
    1/rgamma(1, (length(EM)+1)/2 , (sum(EM^2)+1)/2 )
  }
