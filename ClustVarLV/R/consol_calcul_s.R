consol_calcul_s <-
function(method,X,EXTr,Xr,EXTu,Xu,ind,max.iter=20, eps = 0.01,rlevel)
#
{
  n<-nrow(X)
  p<-ncol(X)
  Xk<-as.matrix(X[,ind])
  pk<-length(ind)
  
  
  dur <- function(a,para)
  {
    b = a
    b[which((abs(a) - para)<0)] = 0
    b
  }
  soft <- function(a,para)
  {    
    #b <- sort(abs(a))?
    b <- abs(a) - para
    b <- (b + abs(b))/2 
    b <- sign(a) * b
    b 
  }
  normvec <- function(a){     
    b = sqrt(sum(a^2))
    if (b==0) b = 1
    b
  }
  trimm <- function(a,para)
  {
    b = a
    b[which((a - para)<0)] = 0
    b
  }
  
  #####################################################################################
  
  if (method==1){
    if ((EXTr==0)&(EXTu==0)) { 
      svd = svd(Xk)
      comp = svd$u[,1]
      CLVcomp = comp * (svd$d[1])
      critere = (svd$d[1])^2/(n-1)
      jj<-which.max(abs(cor(CLVcomp,Xk)))  # modification of the sign of CLVcomp so that it is positively correlated with the closest variable  
      if (cor(CLVcomp,Xk[,jj])<0) CLVcomp=(-1)*CLVcomp
      res<-list(comp=CLVcomp,critere=critere)
    }
#     if((EXTr==1)&(EXTu==0)) {
#       svd = svd(t(Xr)%*%Xk)
#       a = svd$u[,1]
#       comp<- Xr%*% a   
#       CLVcomp = comp
#       critere = (svd$d[1])^2/(n-1)
#       jj<-which.max(cor(CLVcomp,Xk))  # modification of the sign of CLVcomp so that it is positively correlated with the closest variable  
#       if (cor(CLVcomp,Xk[,jj])<0) CLVcomp=(-1)*CLVcomp
#       res<-list(comp=CLVcomp,a=a,critere=critere)     
#     }
#     if ((EXTr==0)&(EXTu==1)) {
#       P<-Xk %*% Xu[ind,]
#       if(sum(P^2)==0) stop("error in P")
#       B<-t(Xk)%*% P
#       vp = eigen(t(B) %*% B)
#       alpha2<-eigen(t(P)%*%P)$values[1] 
#       crit<- vp$values[1]/((n-1)*alpha2)      
#       u<-vp$vectors[,1]
#       comp<-P%*%u /sqrt(alpha2)
#       CLVcomp = comp
#       jj<-which.max(cor(CLVcomp,P))  # modification of the sign of CLVcomp so that it is positively correlated with the closest variable  
#       if (cor(CLVcomp,P[,jj])<0) CLVcomp=(-1)*CLVcomp
#       res<-list(comp=CLVcomp,u=u,critere=crit)
#     }
#     
    
      Ck=comp     
      cova = cov(Xk,Ck)  
     
      para = rlevel*sqrt(diag(var(Xk))/(n-1))
      #beta = dur(cova,para)
      beta = soft(cova,para)
            
      temp <- beta   # in order to check convergence
      temp <- temp/normvec(temp)
      k <- 0
      diff <- 1
      while ((k < max.iter) & (diff > eps)) {
        k <- k + 1 
        
        if ((EXTr==0)&(EXTu==0)){
           alpha <- t(Xk)%*%Xk%*%beta/normvec(t(Xk)%*%Xk%*%beta)
           Ck = Xk%*%alpha/normvec(Xk%*%alpha)
            
        }
#         if ((EXTr==1)&(EXTu==0)){
#           Xrk = t(Xr)%*%Xk
#           alpha <-t(Xrk)%*%Xrk%*%beta/normvec(t(Xrk)%*%Xrk%*%beta)
#           Ck = Xr%*%Xrk%*%alpha/normvec(Xr%*%Xrk%*%alpha)
#         }
#         if ((EXTr==0)&(EXTu==1)){
#           P<-Xk %*% Xu[ind,]
#           if(sum(P^2)==0) stop("error in P")
#           B<-t(Xk)%*% P
#           alpha <- t(B)%*%B%*%beta/normvec(t(B)%*%B%*%beta)
#           Ck = P%*%B%*%alpha/normvec(P%*%B%*%alpha)
#         }
                
        cova = cov(Xk,Ck)                                                
#       beta = dur(cova,para)
        beta = soft(cova,para)

        beta2 = beta/normvec(beta)
        diff <- mean(abs(beta2 - temp))
        temp <- beta2
        
      }
      
      
      beta = beta/normvec(beta)
      colnames(beta) = "loading"   
      res$loading = beta
 
      res$comp = Xk%*%beta     

       if ((EXTr==0)&(EXTu==0)) {
          compnorm=res$comp/normvec(res$comp)
          res$critere<-(n-1)*sum(cov(Xk,compnorm)^2)
        }
      
    
    
#        if ((EXTr==1)&(EXTu==0)) res$critere = sum((Xrk%*%beta)^2)/(n-1)   
        #res$critere = (svd(t(Xrk)%*%Xrk%*%beta)$d[1])^2/(n-1)
#        if ((EXTr==1)&(EXTu==1)) res$critere = sum((B%*%beta)^2)/((n-1)*(svd(P)$d[1])^2)   
        #res$critere = (svd(t(B)%*%B%*%beta)$d[1])^2/(n-1)/((svd(P)$d[1])^2)

      
  }
  ###end of method=1 #########################################################################################
  
  
  
  
  
  #############################################################################################################
  if (method==2){
    if ((EXTu==0)& (EXTr==0)) {
      comp <- Xk %*% matrix(1,pk,1) /pk                 
      #critere<-pk*var(comp)                             # version RSA
      #res<-list(comp=comp,critere=critere)              # version RSA       
      compnorm<-comp/as.numeric(sqrt(t(comp)%*%comp))    # version CommStat  ck normalized
      critere<-pk*apply(comp,2,sd)                       # version CommStat
      alpha = matrix(rep(1/pk,pk),pk,1)
      rownames(alpha)= colnames(Xk)
      res<-list(comp=comp,critere=critere)               # version CommStat
    }    
#     if ((EXTu==0)& (EXTr==1)){                 
#       aa = t(Xr)%*% Xk %*% matrix(1,pk,1) /pk
#       a<-aa/as.numeric(sqrt(t(aa)%*%aa))
#       comp<-(Xr%*%a)
#       critere<-pk*sqrt(t(aa)%*%aa)/(n-1)
#       res<-list(comp=comp,a=a,critere=critere)
#     }
#     if ((EXTu==1)& (EXTr==0)) {
#       Xugroupe<-Xu[ind,]
#       P=Xk%*%Xugroupe
#       if(sum(P^2)==0) stop("error in P")
#       alpha2<-sum(diag(t(P)%*%P))
#       uu = t(P)%*% Xk %*% matrix(1,pk,1) /pk
#       u<-uu/as.numeric(sqrt(t(uu)%*%uu))
#       comp<-(P%*%u)/sqrt(alpha2)
#       critere<-pk*sqrt(t(uu)%*%uu)/sqrt(n-1)
#       critere<-critere/sqrt(alpha2)
#       res<-list(comp=comp,u=u,critere=critere)
#  
#     }   
    
      
  
      
       if ((EXTu==0)& (EXTr==0))  {
         Ck = compnorm
       } else {
         # message. procedure stopped
       }
    
    
      cova = cov(Xk,Ck)
                                   
      para = rlevel*sqrt(diag(var(Xk))/(n-1))
    # beta = dur(cova,para)
      beta = trimm(cova,para)
    
                                                                                   
      temp = beta
      
      k <- 0
      diff <- 1
      while ((k < max.iter) & (diff > eps)) {
        k <- k + 1
        
        if ((EXTr==0)&(EXTu==0)){ 
          # beta includes 0 or 1
          beta[which(beta!=0)]=1
          b = sum(beta)
          if (b==0) b = 1
          alpha <- beta/b
          Cmean = Xk%*%alpha
          Ck = Cmean/normvec(Cmean)  
                                   
        }
#         if ((EXTr==1)&(EXTu==0)){
#           Xrk = t(Xr)%*%Xk
#           alpha <- beta/normvec(beta)
#           aa = Xrk%*%alpha
#           a  = aa/as.numeric(sqrt(t(aa)%*%aa))
#           Cmean = (Xr%*%a)
#           Ck = Cmean/as.numeric(sqrt(t(Cmean)%*%Cmean))       
#         }
#         if ((EXTr==0)&(EXTu==1)){
#           Xugroupe<-Xu[ind,]
#           P=Xk%*%Xugroupe
#           if(sum(P^2)==0) stop("error in P")
#           
#           alpha <- beta/sum(beta)
#           alpha2<-sum(diag(t(P)%*%P))
#           uu = t(P)%*% Xk%*%alpha
#           u = uu/as.numeric(sqrt(t(uu)%*%uu))   
#           Cmean = (P%*%u)/sqrt(alpha2)
#           Ck = Cmean/as.numeric(sqrt(t(Cmean)%*%Cmean))    
#         }
#  
       
        cova = cov(Xk,Ck)
        # beta = dur(cova,para)
        beta = trimm(cova,para)
        if ((EXTr==0)&(EXTu==0))  beta[which(beta!=0)]=1

        beta2 = beta
        diff <- mean(abs(beta2 - temp))
        temp <- beta2
      }
       
      colnames(beta) = "loading"    
      res$comp = Cmean
   #   res$loading = beta
      res$loading =alpha
    
    
        if ((EXTr==0)&(EXTu==0)) {
          compnorm=res$comp/normvec(res$comp)
          res$critere<-sqrt(n-1)*sum(cov(Xk,compnorm))
        }
#        if ((EXTr==1)&(EXTu==0)) res$critere<-length(which(beta!=0))*sd(Cmean) 
#        if ((EXTr==0)&(EXTu==1)) res$critere<-length(which(beta!=0))*sd(uu)/(sqrt(n-1)*sqrt(alpha2))
      
  
  }
  
  return(res)
}
