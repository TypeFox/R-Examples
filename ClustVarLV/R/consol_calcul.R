consol_calcul <-
function(method,X,EXTr,Xr,EXTu,Xu,ind)
#
{
  n<-nrow(X)
  p<-ncol(X)
  Xk<-as.matrix(X[,ind])
  pk<-length(ind)
  
  if (method==1){
    if ((EXTr==0)&(EXTu==0)) { 
       if ( dim(Xk)[1] > dim(Xk)[2] ) {
         vp=eigen(t(Xk)%*%Xk)
#         vp=svd(t(Xk)%*%Xk,nu=1,nv=1)
         comp<-Xk%*%vp$vectors[,1]%*%(vp$values[1])^(-1/2)
#          comp<-Xk%*%vp$u%*%(vp$d)^(-1/2)
         CLVcomp<-Xk%*%vp$vectors[,1]
#          CLVcomp<-Xk%*%vp$u
       }  else {
         vp <- eigen(Xk %*% t(Xk))
#          vp=svd(Xk %*% t(Xk),nu=1,nv=1)
         comp<- vp$vectors[,1]
#          comp<- vp$u
         CLVcomp<- vp$vectors[,1]*sqrt(vp$values[1])
#          CLVcomp<- vp$u*sqrt(vp$d)
       }

      jj<-which.max(cor(CLVcomp,Xk))  # modification of the sign of CLVcomp so that it is positively correlated with the closest variable  
      if (cor(CLVcomp,Xk[,jj])<0) CLVcomp=(-1)*CLVcomp
      critere<-  vp$values[1] /(n-1)
      res<-list(comp=CLVcomp,critere=critere)
    }
    if((EXTr==1)&(EXTu==0)) {
      if ( dim(Xr)[2] > dim(Xk)[2]   ) {  
        vp <- eigen( t(t(Xr)%*%Xk)%*%t(Xr)%*%Xk)
        a<-t(Xr)%*%Xk%*%vp$vectors[,1]%*%(vp$values[1])^(-1/2)
       } else {
        vp <- eigen(t(Xr)%*%Xk %*% t(t(Xr)%*%Xk))
        a<- vp$vectors[,1] 
    }  
                     
      CLVcomp<- Xr%*% a        
      jj<-which.max(cor(CLVcomp,Xk))  # modification of the sign of CLVcomp so that it is positively correlated with the closest variable  
      if (cor(CLVcomp,Xk[,jj])<0) CLVcomp=(-1)*CLVcomp
      critere<-  vp$values[1] /(n-1)
      res<-list(comp=CLVcomp,a=a,critere=critere)
    }
    if ((EXTr==0)&(EXTu==1)) {
      P<-Xk %*% Xu[ind,]
      if(sum(P^2)==0) stop("error in P")
      B<-t(Xk)%*% P
      vp = eigen(t(B) %*% B)
      alpha2<-eigen(t(P)%*%P)$values[1] 
      crit<- vp$values[1]/((n-1)*alpha2)      
      u<-vp$vectors[,1]
      CLVcomp<-P%*%u /sqrt(alpha2)
      jj<-which.max(cor(CLVcomp,P))  # modification of the sign of CLVcomp so that it is positively correlated with the closest variable  
      if (cor(CLVcomp,P[,jj])<0) CLVcomp=(-1)*CLVcomp
      res<-list(comp=CLVcomp,u=u,critere=crit)
    }
  }

  if (method==2){
    if ((EXTu==0)& (EXTr==0)) {
      comp <- Xk %*% matrix(1,pk,1) /pk                 
      #critere<-pk*var(comp)                             # version RSA
      #res<-list(comp=comp,critere=critere)              # version RSA       
      compnorm<-comp/as.numeric(sqrt(t(comp)%*%comp))  # version CommStat  ck normalized
      critere<-pk*apply(comp,2,sd)                             # version CommStat
      res<-list(comp=comp,critere=critere)         # version CommStat
    }    
    if ((EXTu==0)& (EXTr==1)){                 
      aa = t(Xr)%*% Xk %*% matrix(1,pk,1) /pk
      a<-aa/as.numeric(sqrt(t(aa)%*%aa))
      comp<-(Xr%*%a)
      critere<-pk*sqrt(t(aa)%*%aa)/(n-1)
      res<-list(comp=comp,a=a,critere=critere)
    }
    if ((EXTu==1)& (EXTr==0)) {
      Xugroupe<-Xu[ind,]
      P=Xk%*%Xugroupe
      if(sum(P^2)==0) stop("error in P")
      alpha2<-sum(diag(t(P)%*%P))
      uu = t(P)%*% Xk %*% matrix(1,pk,1) /pk
      u<-uu/as.numeric(sqrt(t(uu)%*%uu))
      comp<-(P%*%u)/sqrt(alpha2)
      critere<-pk*sqrt(t(uu)%*%uu)/sqrt(n-1)
      critere<-critere/sqrt(alpha2)
      res<-list(comp=comp,u=u,critere=critere)
 
    }   
  }
  
  return(res)
}
