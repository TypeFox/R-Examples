# @title Regression coefficients conversion
# @description converts the regression coefficients associated with the components 
# into coefficients associated with the initial regressors
# @export 
# @param Xcr matrix of the standardized regressors
# @param centerx vector of the regressors' means
# @param invsqrtm metric matrix M^(-1/2)
# @param gamma (generalized) regression coefficients associated with the components 
# @param u matrix of the component coordinates on Xcr (component = Xcr * u) 
# @param comp matrix of the components
# @return beta matrix containing the regression coefficients associated with the regressors X
f2x <- function(Xcr,centerx,invsqrtm,gamma,u,comp) {
  n <- dim(Xcr)[1]
  p<-dim(Xcr)[2]
  q<- dim(gamma)[2]
  K <- ncol(comp)  
  gamma <- as.matrix(gamma)
  if (K>1) {
    pi <- diag(p)
    V<- matrix(u[,1],p,1)
    for (k in 2:K) {      
      pi<-pi-( tcrossprod(V[,k-1,drop=F],comp[,k-1,drop=F])%*%Xcr%*%pi)/sum(comp[,k-1]^2)
      V <- cbind(V,pi%*%u[,k])
    }
    beta <- invsqrtm%*%V%*%gamma[-1,,drop=F]##
    alpha<-gamma[1,]-centerx%*%beta
    beta<-rbind(alpha,beta)
  } else {
    beta<-invsqrtm%*%u%*%gamma[-1,,drop=FALSE]##
    alpha<-gamma[1,]-centerx%*%beta
    beta<-rbind(alpha,beta)
  }
  return(beta)
}
