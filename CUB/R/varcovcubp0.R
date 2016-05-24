# @title Variance-covariance matrix of CUB model with covariates for the uncertainty parameter
# @description Compute the variance-covariance matrix of parameter estimates of a CUB model with 
# covariates for the uncertainty component.
# @aliases varcovcubp0
# @usage varcovcubp0(m, ordinal, Y, bet, csi)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of covariates for explaining the uncertainty parameter
# @param bet Vector of parameters for the uncertainty component, whose length equals NCOL(Y)+1 
# to include an intercept term (first entry)
# @param csi Feeling parameter
# @details The function checks if the variance-covariance matrix is positive-definite: if not, 
# it returns a warning message and produces a matrix with NA entries.
# @seealso \code{\link{probcubpq}}, \code{\link{cubpq}}
# @references
# Piccolo D. (2006), Observed Information Matrix for CUB Models, \emph{Quaderni di Statistica}, \bold{8}, 33--78
#' @keywords internal

varcovcubp0 <-
function(m,ordinal,Y,bet,csi){
  vvi<-(m-ordinal)/csi-(ordinal-1)/(1-csi)
  ui<-(m-ordinal)/(csi^2)+(ordinal-1)/((1-csi)^2)
  qi<-1/(m*probcubp0(m,ordinal,Y,bet,csi))
  ei<-logis(Y,bet)
  qistar<-1-(1-ei)*qi
  eitilde<-ei*(1-ei)
  qitilde<-qistar*(1-qistar)
  ff<-eitilde-qitilde
  g10<-vvi*qitilde
  YY<-cbind(1,Y)                           
  i11<-t(YY)%*%(Hadprod(YY,ff))            ###ALTERNATIVE  YY*ff does not work
  i12<- -t(YY)%*%(g10) 
  i22<-sum(ui*qistar-(vvi^2)*qitilde)
  # Information matrix 
  ## matinf=rbind(cbind(i11,i12),cbind(t(i12),i22))
  # Var-covar matrix
  
  nparam<- NCOL(Y) + 2
  matinf<-matrix(NA,nrow=nparam,ncol=nparam)
  for (i in 1:(nparam-1)){
    matinf[i,]<-t(c(i11[i,],i12[i]))
  }
  matinf[nparam,]<-t(c(t(i12),i22))
  
  
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
