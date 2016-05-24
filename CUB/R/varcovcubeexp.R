# @title Variance-covariance matrix for CUBE models based on the expected information matrix
# @description Compute the variance-covariance matrix of parameter estimates as the inverse of
#  the expected information matrix for a CUBE model without covariates.
# @aliases varcovcubeexp
# @usage varcovcubeexp(m, pai, csi, phi, n)
# @param m Number of ordinal categories
# @param pai Uncertainty parameter
# @param csi Feeling parameter
# @param phi Overdispersion parameter
# @param n Number of observations
# @details The function checks if the variance-covariance matrix is positive-definite: if not, 
# it returns a warning message and produces a matrix with NA entries.
# @seealso \code{\link{cube000}}, \code{\link{varcovcubeobs}}
# @references
# Iannario, M. (2014). Modelling Uncertainty and Overdispersion in Ordinal Data, 
# \emph{Communications in Statistics - Theory and Methods}, \bold{43}, 771--786
#' @keywords internal


varcovcubeexp <-
function(m,pai,csi,phi,n){
  pr<-probcube(m,pai,csi,phi);
  sum1<-sum2<-sum3<-sum4<-rep(NA,m);
  for(jr in 1:m){
    seq1<-1/((1-csi)+phi*((1:jr)-1));
    seq2<-1/((csi)+phi*((1:(m-jr+1))-1));
    seq3<-((1:jr)-1)/((1-csi)+phi*((1:jr)-1));
    seq4<-((1:(m-jr+1))-1)/((csi)+phi*((1:(m-jr+1))-1));
    sum1[jr]<-sum(seq1)-seq1[jr];
    sum2[jr]<-sum(seq2)-seq2[m-jr+1];
    sum3[jr]<-sum(seq3)-seq3[jr];
    sum4[jr]<-sum(seq4)-seq4[m-jr+1];
  }
  sum5<-sum(((1:m)-1)/(1+phi*(1:m)));
  
  cr<-pr-(1-pai)/m;
  derpai<-(pr-1/m)/pai;
  dercsi<-cr*(sum2-sum1);
  derphi<-cr*(sum3+sum4-sum5);
  
  infpaipai<-sum((derpai^2)/pr);
  infpaicsi<-sum((derpai*dercsi)/pr);
  infpaiphi<-sum((derpai*derphi)/pr);
  infcsicsi<-sum((dercsi^2)/pr);
  infcsiphi<-sum((dercsi*derphi)/pr);
  infphiphi<-sum((derphi^2)/pr);
  ### Information matrix 
  inform<-matrix(NA,nrow=3,ncol=3);
  inform[1,1]<-infpaipai;   inform[1,2]<-infpaicsi;   inform[1,3]<-infpaiphi;
  inform[2,1]<-infpaicsi;   inform[2,2]<-infcsicsi;   inform[2,3]<-infcsiphi;
  inform[3,1]<-infpaiphi;   inform[3,2]<-infcsiphi;   inform[3,3]<-infphiphi;
  ### Var-covar matrix 
  
  if(any(is.na(inform))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=3,ncol=3)
  } else {
    if(det(inform)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=3,ncol=3)
    } else {
      varmat<-solve(inform)/n
    }
  }
  
  return(varmat)
}
