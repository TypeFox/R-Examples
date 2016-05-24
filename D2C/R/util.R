
rankrho<-function(X,Y,nmax=5,regr=FALSE,first=NULL){
  ## mutual information ranking
  ## 17/10/11
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y)

  if (var(Y)<0.01)
    return(1:nmax)
  X<-scale(X)
  
  Iy<-numeric(n)
  if (!regr){
    Iy<-cor2I2(corXY(X,Y))
  } else {
    for (i in 1:n)
      Iy[i]<-abs(regrlin(X[,i],Y)$beta.hat[2])
  }
  
  if (m>1)
    Iy<-apply(Iy,1,mean)
  
  
  return(sort(c(Iy), decreasing=TRUE, index.return=TRUE)$ix[1:nmax])
  
  
}




quantization<-function(x,nbin=1){
  if (nbin==1)
    return(as.numeric(cut(x, breaks = c(min(x)-1,median(x),max(x)+1))))
  else
    return(as.numeric(cut(x, breaks = c(min(x)-1,quantile(x,c(0.25,0.5,0.75)),max(x)+1))))

}

H_sigmoid <- function(n=2)
{
  a = runif(n+1,min = -5,max = 5)
  f <- function(x)
  {
    X  = x^(0:n)
    return ( 1/(1+exp(sum(X * a))))
  }
  return(Vectorize(f))
}
H_Rn <- function(n){
  a = runif(n+1,min = -1,max = 1)
  f <- function(x)
  {
    X  = x^(0:n)
    return ( sum(X * a))
  }
  return(Vectorize(f))
}


pcor1<-function(x,y,z){
  ## partial correlation cor(x,y|z)

  
  if (is.numeric(z)){
    rho.xy<-cor(x,y,"pairwise.complete.obs")
    rho.xz<-cor(x,z,"pairwise.complete.obs")
    rho.yz<-cor(z,y,"pairwise.complete.obs")
    if (is.na(rho.xz+rho.yz+rho.xy))
      return(0)
    if (rho.xz==1 | rho.yz==1)
      return(0)
    rho<-(rho.xy-rho.xz*rho.yz)/(sqrt(1-min(rho.xz^2,0.99))*sqrt(1-min(rho.yz^2,0.99)))
    return(rho)
  } else {
    stop("z should be numeric")
    
  }
  
  
  
}


corDC<-function(X,Y){
  ## correlation continuous matrix and discrete vector
  ## NB: the notion of sign has no meaning in this case. Mean of absolute values is taken
  ## 14/11/2011
  
  if (!is.factor(Y))
    stop("This is not the right function. Y is not a factor !!")
  
  N<-NROW(X)
  L<-levels(Y)
  
  if( length(L)==2)
    lL<-1
  else
    lL<-length(L)
  
  cxy<-NULL
  for (i in 1:lL){
    yy<-numeric(N)
    ind1<-which(Y==L[i])
    ind2<-setdiff(1:N,ind1)
    yy[ind1]<-1
    cxy<-cbind(cxy,abs(cor(X,yy)))
  }
  
  apply(cxy,1,mean)
}


Icond<-function(x,y=NULL,z,lambda=0){
  ## conditional  information cor(x,y|z)

  
  
  ## numeric z
  if (is.numeric(z)){
    if (is.vector(x))
      return(cor2I2(pcor1(x,y,z)))
    X<-x
    n<-NCOL(X)
    Ic<-array(0,c(n,n))
    for (i in 1:(n-1))
      for (j in (i+1):n){
        Ic[i,j]<-Icond(X[,i],X[,j],z)
        Ic[j,i]<-Ic[i,j]
      }
    return(Ic)
    
  }
  ## factor z and vectors x and y 
  if (! is.null(y)){
    L<-levels(z)
    lL<-length(L)
    w<-numeric(lL)
    for (i in 1:lL)
      w[i]<-length(which(z==L[i]))
    w<-w/sum(w)
    
    Ic<-NULL
    for (i in 1:lL){
      ind1<-which(z==L[i])
      Ic<-c(Ic,cor2I2(cor(x[ind1],y[ind1])))
    }
    
    return(as.numeric(w*Ic))
  }
  
  ## factor z and matrix x
  X<-x
  n<-NCOL(X)
  L<-levels(z)
  lL<-length(L)
  w<-numeric(lL)
  for (i in 1:lL)
    w[i]<-length(which(z==L[i]))
  w<-w/sum(w)
  
  Ic<-array(0,c(n,n))
  W<-0
  for (i in 1:lL){
    ind1<-which(z==L[i])
    
    
    if (length(ind1)>8){
      
      
      Ic<-Ic+w[i]*cor2I2(cor.shrink(X[ind1,],lambda=lambda,verbose=F))
      W<-W+w[i]
    }
  }
  
  
  return(Ic/W)
  
}







  
  
ppears<-function(r.hat,N,S=0){
  n<-length(r.hat)
  p<-numeric(n)
  
  for (i in 1:n){
    z<-abs(0.5*(log(1+r.hat[i])-log(1-r.hat[i])))*sqrt(N[i]-S-3)
    
    p[i]<-pnorm(z,lower.tail=F)
    
    
  }
  p
}

corXY<-function(X,Y){
  ## correlation continuous matrix and continuous/discrete vectormatrix
  
  
  n<-NCOL(X)
  N<-NROW(X)
  m<-NCOL(Y)
  
  cXY<-array(NA,c(n,m))
  
  for (i in 1:m){
    if (m==1)
      YY<-Y
    else
      YY<-Y[,i]
    if (is.numeric(YY)){
      cXY[,i]<-cor(X,YY,use="pairwise.complete.obs")
    } else {
      cXY[,i]<-corDC(X,YY)
    }
  }
  cXY
}


cor2I2<-function(rho){
  rho<-pmin(rho,1-1e-5)
  rho<-pmax(rho,-1+1e-5)
  -1/2*log(1-rho^2)
  
  
}


lazy.pred<- function(X,Y,X.ts,class=FALSE,return.more=FALSE,
                     conPar=3,linPar=5,cmbPar=10){

  n<-NCOL(X)
  N<-NROW(X)
  
  if (class){ ## classification
    l.Y<-levels(Y)
    L<-length(l.Y)
    u<-unique(Y)
    
    if (length(u)==1){
      P<-array(0,c(NROW(X.ts),L))
      colnames(P)<-l.Y
      P[,u]<-1
      out.hat<-factor(rep(as.character(u),length(X.ts)),levels=l.Y)
      return(list(pred=out.hat,prob=P))
    }
    
    if (L==2) {
      
      
      stop("not supported")
    } else {
      algo="lazy"
      
      stop("not supported")
      
    }
  } else { ## regression
    d<-data.frame(cbind(Y,X))
    names(d)[1]<-"Y"
    names(d)[2:(n+1)]<-paste("x",1:n,sep="")
    
    
    
    mod<-lazy(Y~.,d,control=lazy.control(distance="euclidean",
                                         conIdPar=conPar,
                                         linIdPar=linPar,
                                         cmbPar=cmbPar))
    if (is.vector(X.ts) & n>1)
      X.ts<-array(X.ts,c(1,n))
    d.ts<-data.frame(X.ts)
    
    names(d.ts)<-names(d)[2:(n+1)]
    
    if (!return.more){
      ll<- predict(mod,d.ts)
      return(ll$h)
      
    } else {
      ll<- predict(mod,d.ts,S.out=TRUE,k.out=FALSE)
      return(ll)
    }
  }
  
}



#' mIMR (minimum Interaction max Relevance) filter 
#' @param X :  input matrix 
#' @param Y : output vector
#' @param nmax : number of returned features
#' @param init : if TRUE it makes a search in the space of pairs of features to initialize the ranking, otherwise the first ranked feature is the one with the highest mutual information with the output
#' @param lambda : weight \eqn{0 \le \lambda \le 1} of the interaction term
#' @param spouse.removal : TRUE OR FALSE. if TRUE it removes the spouses before ranking
#' @param caus :   if \code{caus =1} it prioritizes causes otherwise (\code{caus=-1}) it prioritizes effects
#' @return ranked vector of \code{nmax} indices of features
#' @description Filter  based on information theory which aims to prioritise direct causal relationships in feature selection problems where the ratio between the number of features and the number of samples is high. The approach is based on the notion of interaction which is informative about the relevance of an input subset as well as its causal relationship with the target. 
#' @examples
#'  set.seed(0)
#' N<-500
#' n<-5
#' X<-array(rnorm(N*n),c(N,n))
#' Y<-X[,1]-3*X[,3]+4*X[,2]+rnorm(N,sd=0.5)
#' Z1<-Y+rnorm(N,sd=0.5)
#' ## effect 1
#' Z2<-2*Y+rnorm(N,sd=0.5)
#' ## effect 2
#' most.probable.causes<-mimr(cbind(X,Z1,Z2),Y,nmax=3,init=TRUE,spouse=FALSE,lambda=1)
#' ## causes are in the first three columns of the feature dataset 
#' most.probable.effects<-mimr(cbind(X,Z1,Z2),Y,nmax=3,init=TRUE,spouse=FALSE,lambda=1,caus=-1)
#' ## effects are in the last two columns of the feature dataset
#' @references Bontempi G., Meyer P.E. (2010) Causal filter selection in microarray data. ICML10 
#' @export
mimr<-function(X,Y,nmax=5,
               init=FALSE,lambda=0.5,
               spouse.removal=TRUE,
               caus=1){
  if (var(Y)<0.01)
    return(1:nmax)
  NMAX<-nmax
  m<-NCOL(Y) # number of outputs
  n<-NCOL(X)
  orign<-n
  N<-NROW(X)
  H<-apply(X,2,var)
  HY<-var(Y)
  CY<-corXY(X,Y)
  Iy<-cor2I2(CY)  
  subset<-1:n
  pv.rho<-ppears(c(CY),N+numeric(n))
  if (spouse.removal){
    pv<-ppears(c(CY),N+numeric(n))
    s<-sort(pv,decreasing=TRUE,index.return=TRUE)  
    hw<-min(n-nmax,length(which(s$x>0.05)))
    spouse<-s$ix[1:hw]
    subset<-setdiff(1:n,s$ix[1:hw])
    X<-X[,subset]
    Iy<-Iy[subset]
    n<-NCOL(X)
  }
  
  
  CCx<-cor(X) 
  Ix<-cor2I2(CCx)
  ## mutual information
  Ixx<-Icond(X,z=Y,lambda=0.02)
  ## conditional information
  Inter<-array(NA,c(n,n))
 
  if (init){    
    max.kj<--Inf
    for (kk in 1:(n-1)){
      for (jj in (kk+1):n){
        Inter[kk,jj]<- (1-lambda)*(Iy[kk]+Iy[jj])+caus*lambda*(Ixx[kk,jj]-Ix[kk,jj])
        Inter[jj,kk]<-Inter[kk,jj]
        if (Inter[kk,jj]>max.kj){
          max.kj<-Inter[kk,jj]
          subs<-c(kk,jj)
        }
      }
    }
  } else {
    subs<-which.max(Iy)
  }
  
  if (nmax>length(subs)){
    last.subs<-0
    for (j in length(subs):min(n-1,NMAX-1)){
      mrmr<-numeric(n)-Inf
      if (length(subs)<(n-1)){
        if (length(subs)>1){
          mrmr[-subs]<- (1-lambda)*Iy[-subs]+caus*lambda*apply(-Ix[subs,-subs]+Ixx[subs,-subs],2,mean)
        } else {
          mrmr[-subs]<- (1-lambda)*Iy[-subs]+caus*lambda*(-Ix[subs,-subs]+Ixx[subs,-subs])
        }
      } else {
        mrmr[-subs]<-Inf
      }
      s<-which.max(mrmr)
      subs<-c(subs,s)  
    }
    
    
    
    
  }
  
  ra<-subset[subs]
  
  if (nmax>length(ra))
    ra<-c(ra,setdiff(1:orign,ra))
  
  ra
  
}


assoc <-function(x,y){
  c(abs(cor(x,y)),cor.test(x,y)$p.value)
 
}


#' Balanced Error Rate
#' @param Ytrue :  binary numeric vector (made of 0 or 1) of real classes 
#' @param Yhat : binary numeric vector (made of 0 or 1) of predicted classes
#' @description The balanced error rate is the average of the errors on each class: BER = 0.5*(FP/(TN+FP) + FN/(FN+TP)).
#' @return Balanced Error Rate \eqn{0 \le } BER \eqn{ \le 1}
#' @export
BER<-function(Ytrue,Yhat){

  if (!(is.numeric(Ytrue) & is.numeric(Yhat)))
    stop("BER accepts only numeric values")
  TN<-length(which(Yhat==0 & Ytrue==0)) 
  FN<-length(which(Yhat==0 & Ytrue==1))
  TP<-length(which(Yhat==1 & Ytrue==1))
  FP<-length(which(Yhat==1 & Ytrue==0))

  
  b1<-FP/(TN+FP)
  b2<-FN/(FN+TP)
  if (is.na(b1))
    b1<-0
  if (is.na(b2))
    b2<-0
  return(0.5*(b1+b2))
  


}


