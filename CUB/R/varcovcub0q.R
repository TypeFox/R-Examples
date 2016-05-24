# @title Variance-covariance matrix of CUB models with covariates for the feeling component
# @description Compute the variance-covariance matrix of parameter estimates of a CUB model
#  with covariates for the feeling component.
# @aliases varcovcub0q
# @usage varcovcub0q(m, ordinal, W, pai, gama)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param W Matrix of covariates for explaining the feeling component
# @param pai Uncertainty parameter
# @param gama Vector of parameters for the feeling component, whose length is 
# NCOL(W)+1 to include an intercept term in the model (first entry of gama)
# @export varcovcub0q
# @details The function checks if the variance-covariance matrix is positive-definite: if not, 
# it returns a warning message and produces a matrix with NA entries.
# @seealso \code{\link{probcub0q}}, \code{\link{loglikcub0q}}, \code{\link{cub0q}}  
# @references
# Piccolo D.(2006), Observed Information Matrix for MUB Models. \emph{Quaderni di Statistica},
#  \bold{8}, 33--78,
# @examples
# data(univer)
# m<-7
# ordinal<-univer[,9]
# pai<-0.86
# gama<-c(-1.94, -0.17)
# W<-univer[,4]           
# varmat<-varcovcub0q(m, ordinal, W, pai, gama)
#' @keywords internal


varcovcub0q <-
function(m,ordinal,W,pai,gama){
  qi<-1/(m*probcub0q(m,ordinal,W,pai,gama))
  qistar<-1-(1-pai)*qi
  qitilde<-qistar*(1-qistar)
  fi<-logis(W,gama)
  fitilde<-fi*(1-fi)
  ai<-(ordinal-1)-(m-1)*(1-fi)
  g01<-(ai*qi*qistar)/pai 
  hh<-(m-1)*qistar*fitilde-(ai^2)*qitilde
  WW<-cbind(1,W)                         
  i11<-sum((1-qi)^2)/(pai^2) 
  i12<-t(g01)%*%WW 
  i22<-t(WW)%*%(Hadprod(WW,hh))                    ### i22=t(WW)%*%(WW*hh) does not work;  
  
  ### Information matrix 
  ##matinf=rbind(cbind(i11,i12),cbind(t(i12),i22))
  nparam<-NCOL(W)+2
  
  matinf<-matrix(NA,nrow=nparam,ncol=nparam)
  matinf[1,]<-t(c(i11,i12))
  for (i in 2:(nparam)){
    matinf[i,]<-t(c(i12[i-1],i22[i-1,]))
  }
  
  if(any(is.na(matinf))==TRUE){
    warning("ATTENTION: NAs produced")
    varmat<-matrix(NA,nrow=nparam,ncol=nparam)
  } else {
    if(det(matinf)<=0){  
      warning("ATTENTION: Variance-covariance matrix NOT positive definite")
      varmat<-matrix(NA,nrow=nparam,ncol=nparam)
    } else {
      varmat<-solve(matinf)
    }
  }
  
  return(varmat)
}
