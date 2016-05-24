cal_df_gamma1 <-
function(yXb, poly.index, h, nbins, plmm){
  bing<-binning(y=yXb, x=as.vector(plmm$T_mat), nbins=nbins)
  x_til<-bing$x
  M<-length(x_til)
  C<-matrix(rep(bing$x.freq, M), ncol=M, byrow=T) 
## = rep(1, M)%x%t(bing$x.freq)
  X_til<-matrix(rep(x_til, rep(M, M)), byrow=T, ncol=M) - matrix(rep(x_til, rep(M, M)), byrow=F, ncol=M)
## = x_til%x%t(rep(1, M))-rep(1, M)%x%t(x_til)
## symmetric
## 1.row of X = elements of x-x.g1 
## 2.row of X = elements of x-x.g2 etc

  W_til<-exp(-0.5*(X_til/h)^2)
  WC<-W_til*C
## 1.row of WC = elements of W.g1%*%C 
## 2.row of WC = elements of W.g2%*%C  etc    
  
### tr(S) ###     
  if(poly.index==1){### Local Linear
    s0<-WC%*%rep(1,M)## vector {s0.l} l=1~M
    s1<-(WC*X_til)%*%rep(1,M)## vector {s1.l} l=1~M
    s2<-(WC*X_til^2)%*%rep(1,M)## vector {s2.l} l=1~M

    s0s2_s1s1<-s0*s2-s1^2## vector
## determinant for (XWCX)^(-1)
    S.ll<-s2/s0s2_s1s1## diag(S_til) vector
    tr_S<-sum(S.ll*bing$x.freq)
  }else if (poly.index==0){
    S.ll<-1/(WC%*%rep(1,M))## XWCX = scalar
    tr_S<-sum(S.ll*bing$x.freq)
  }
  
### tr(SS') ###
  if(poly.index==1){### Local Linear
### XWWCX
    WWC<-W_til*WC
## 1.row of WWC = elements of W.g1^2%*%C 
## 2.row of WWC = elements of W.g2^2%*%C etc

    z0<-WWC%*%rep(1,M)## vector {z0.l} l=1~M
    z1<-(WWC*X_til)%*%rep(1,M)## vector {z1.l} l=1~M
    z2<-(WWC*X_til^2)%*%rep(1,M)## vector {z2.l} l=1~M

### e(XWCX)XWWCX(XWCX)e 
    a.l<-s2/s0s2_s1s1
    b.l<--s1/s0s2_s1s1
    SS.ll<-a.l^2*z0+2*a.l*b.l*z1+b.l^2*z2
    tr_SS<-sum(SS.ll*bing$x.freq)
  }else if (poly.index==0){### NW
    WWC<-W_til*WC
    SS.ll<-S.ll^2*(WWC%*%rep(1,M))
    tr_SS<-sum(SS.ll*bing$x.freq)
  }     
  return(2*tr_S-tr_SS)  
}
