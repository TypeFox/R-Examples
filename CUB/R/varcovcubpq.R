# @title Variance-covariance matrix of a CUB model with covariates for both uncertainty and feeling
# @description Compute the variance-covariance matrix of parameter estimates of a CUB model with covariates for
#  both the uncertainty and the feeling components.
# @aliases varcovcubpq
# @usage varcovcubpq(m, ordinal, Y, W, bet, gama)
# @param m Number of ordinal categories
# @param ordinal Vector of ordinal responses
# @param Y Matrix of covariates for explaining the uncertainty parameter
# @param W Matrix of covariates for explaining the feeling parameter
# @param bet Vector of parameters for the uncertainty component, with length equal to
# NCOL(Y)+1 to account for an intercept term (first entry)
# @param gama Vector of parameters for the feeling component, with length equal to
# NCOL(W)+1 to account for an intercept term (first entry)
# @details The function checks if the variance-covariance matrix is positive-definite: if not, it returns a warning
#  message and produces a matrix with NA entries.
# @seealso \code{\link{probcubpq}}, \code{\link{cubpq}}
# @references
# Piccolo D. (2006), Observed Information Matrix for CUB Models, \emph{Quaderni di Statistica}, \bold{8}, 33--78
#' @keywords internal

varcovcubpq <-
function(m,ordinal,Y,W,bet,gama){
  qi<-1/(m*probcubpq(m,ordinal,Y,W,bet,gama))
  ei<-logis(Y,bet)
  eitilde<-ei*(1-ei)
  qistar<-1-(1-ei)*qi
  qitilde<-qistar*(1-qistar)    
  fi<-logis(W,gama)
  fitilde<-fi*(1-fi)
  ai<-(ordinal-1)-(m-1)*(1-fi)
  ff<-eitilde-qitilde
  gg<-ai*qitilde
  hh<-(m-1)*qistar*fitilde-(ai^2)*qitilde
  YY<-cbind(1,Y)                      
  WW<-cbind(1,W)                       
  i11<-t(YY)%*%(Hadprod(YY,ff))    ###   i11<-t(YY)%*%(YY*ff); 
  i12<-t(YY)%*%(Hadprod(WW,gg))    ###   i12<-t(YY)%*%(WW*gg); 
  i22<-t(WW)%*%(Hadprod(WW,hh))    ###   i22<-t(WW)%*%; 
  ## matinf<-rbind(cbind(i11,i12),cbind(t(i12),i22))  # Information matrix 
  nparam<-NCOL(Y)+NCOL(W)+2
  
  npai<-NCOL(Y)+1
  ncsi<-NCOL(W)+1
  nparam<-npai+ncsi
  
  matinf<-matrix(NA,nrow=nparam,ncol=nparam)
  for (i in 1:npai){
    matinf[i,]<-t(c(i11[i,],i12[i,]))
  }
  for (i in (npai+1):nparam){
    matinf[i,]<-t(c(t(i12)[i-npai,],i22[i-npai,]))
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
