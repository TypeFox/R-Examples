cal_df_gamma2 <-
function(yXb, poly.index, h, nbins, plmm){
  bing<-binning(y=yXb,x=plmm$T_mat,nbins=nbins)
  x_til<-bing$x
  M<-nrow(x_til)
  C<-matrix(rep(bing$x.freq, M), ncol=M, byrow=T) 
## = rep(1, M)%x%t(bing$x.freq)
  X_tilde1<-matrix(rep(x_til[,1], rep(M, M)), byrow=T, ncol=M) - matrix(rep(x_til[,1], rep(M, M)), byrow=F, ncol=M)
## = x_til%x%t(rep(1, M))-rep(1, M)%x%t(x_til)
## symmetric
## 1.row of X = x.g1-x.gl (l=1~M) 
## 2.row of X = x.g2-x.gl etc
  W_tilde<-exp(-0.5*(X_tilde1/h[1])^2)
  
  X_tilde2<-matrix(rep(x_til[,2], rep(M, M)), byrow=T, ncol=M) - matrix(rep(x_til[,2], rep(M, M)), byrow=F, ncol=M)
  W_tilde<-W_tilde*exp(-0.5*(X_tilde2/h[2])^2)
## Multiplicative Kernel
  
  WC<-W_tilde*C
## 1.row of WC = elements of W.g1%*%C 
## 2.row of WC = elements of W.g2%*%C  etc    
  
### tr(S) ###     
  if(poly.index==1){### Local Linear
    s0<-WC%*%rep(1,M)## vector {s0.l} l=1~M
    s1<-(WC*X_tilde1)%*%rep(1,M)## vector {s1.l} l=1~M
    s2<-(WC*X_tilde2)%*%rep(1,M)
    s11<-(WC*X_tilde1^2)%*%rep(1,M)
    s12<-(WC*X_tilde1*X_tilde2)%*%rep(1,M)
    s22<-(WC*X_tilde2^2)%*%rep(1,M)
    
## determinant for (XWCX)^(-1). this is vector    
    det<-(s0*s11*s22)+2*(s1*s12*s2)-s11*s2^2-s0*s12^2-s1^2*s22
    
    S.ll<-(s11*s22-s12^2)/det## diag(S_til) vector   
    tr_S<-sum(S.ll*bing$x.freq)
  }else if (poly.index==0){### Nadaraya-Watson
    S.ll<-1/(WC%*%rep(1,M))## XWCX = scalar
    tr_S<-sum(S.ll*bing$x.freq)
  }
  
### tr(SS') ###
  if(poly.index==1){### Local Linear
### XWWCX
    WWC<-W_tilde*WC
## 1.row of WWC = elements of W.g1^2%*%C 
## 2.row of WWC = elements of W.g2^2%*%C etc

    z0<-WWC%*%rep(1,M)## vector {z0.l} l=1~M
    z1<-(WWC*X_tilde1)%*%rep(1,M)## vector {s1.l} l=1~M
    z2<-(WWC*X_tilde2)%*%rep(1,M)
    z11<-(WWC*X_tilde1^2)%*%rep(1,M)
    z12<-(WWC*X_tilde1*X_tilde2)%*%rep(1,M)
    z22<-(WWC*X_tilde2^2)%*%rep(1,M)

### e(XWCX)XWWCX(XWCX)e 
    a.l<-(s11*s22-s12^2)/det
    b.l<--(s1*s22-s2*s12)/det
    c.l<-(s1*s12-s2*s11)/det
    
    SS.ll<-a.l^2*z0+a.l*b.l*z1+a.l*c.l*z2+
          a.l*b.l*z1+b.l^2*z11+b.l*c.l*z12+
          a.l*c.l*z2+b.l*c.l*z12+c.l^2*z22        
          
    tr_SS<-sum(SS.ll*bing$x.freq)
  }else if (poly.index==0){### NW
    WWC<-W_tilde*WC
    SS.ll<-S.ll^2*(WWC%*%rep(1,M))
    tr_SS<-sum(SS.ll*bing$x.freq)
  }   
  
  return(2*tr_S-tr_SS)  
}
