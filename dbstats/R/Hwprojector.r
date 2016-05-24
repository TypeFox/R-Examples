  
  
 ######################
 #### Projector Hw ####
 ######################

 ## Desciption: Calculates and returns the matrix of the projector Hw.
 ##                   " Hw <- G %*% Dsqw %*% Fwplus %*% Dsqw ".
 ##         Where Dw=diag(weights) is diagonal matrix whose diagonal entries are 
 ##         the weights weights and Dsqw = sqrt(Dw). The matrix G is the n x n 
 ##         inner-products matrix and finally, Fwplus is the Moore-Penrose 
 ##         pseudo-inverse of Fw.
 ##
 ##         By performing the inverse Fw, its running a singular value 
 ##         decomposition for the Fw, such that Fw=UDeltaU'. It's easy to find 
 ##         Fwplus as U(1/Delta)U'.
 ##
 ##         The key point is that the inverse is not the whole matrix. 
 ##         We build the mentioned Fwplus with eff.rank components which ones are 
 ##         chosen from a certain method.
 ## 

 HwProject<-function(G,Dsqw,rk,epsilon,cvyes=FALSE,n){
   
   Fw <- Dsqw %*% G %*% Dsqw

   # cvyes=true indicate that a defined method, like "OCV", is used. 
   # Find the rank that minimizes the statistic defined by the method used)
   # (for all the effective ranks to explain 99% of the variability.
   
   epsilon_CV<-1.0e-9  

   # Trim to rank ri. makes the singular value descomposition of Fw.
   Fw.svd<-svd(Fw)    
   Fw.u<-Fw.svd$u
   Fw.v<-Fw.svd$v
   Fw.d<-Fw.svd$d
   
   # Singular values in percent
   Fw.dpcent<-Fw.d/sum(Fw.d)   
   
   # acumulative singular values in percent 
   Fw.dnou<-apply(as.matrix(1:length(Fw.d)),1,function(i){
                      sum(Fw.dpcent[1:i])  
      })

   # for eff.rank method (not methods ocv, gcv, aic and bic).
   # Take as many components as eigenvalues (Fw.d) "exceeding" epsilon.
   if (!cvyes){    
    if (rk == 0){
     eff.rank<-
       if (any((1-round(Fw.dnou,10))<=epsilon)) 
          which((1-round(Fw.dnou,10))<=epsilon)[1] 
       else 
          1   # si l'eff.rank==0 --> s'agafa el 1-epsilon% de variabilitat
     
     if (eff.rank==n) eff.rank<-n-1 # si l'eff.rank==n --> s'agafa n-1 components
    }else{
       Fw.r<-
       if (any((1-round(Fw.dnou,10))<=epsilon_CV)) 
          which((1-round(Fw.dnou,10))<=epsilon_CV)[1] 
       else 
          1   # si l'eff.rank==0 --> s'agafa el 1-epsilon% de variabilitat
       
       # if selected rk is to higher --> default method : Fw.r<-sum(Fw.d>epsilon) 
       if (rk > Fw.r){      
        warning(sprintf("  maximum effective rank = %i, lower than the 'eff.rank' used = %i. Will apply the default parameters",Fw.r,rk)) 
       }
       eff.rank <- min(rk,Fw.r) 
    }  

    # use methods ocv, gcv, aic or bic 
   }else{      
     # threshold (effective rank maxim) that the methods computing the HHat matrix 
     # (at maximum, takes as many components such that at least, takes the 
     # 99% variability of the data.
     Fw.r<-
       if (any((1-round(Fw.dnou,10))<=epsilon_CV)) 
          which((1-round(Fw.dnou,10))<=epsilon_CV)[1] 
       else 
          1   # si l'eff.rank==0 --> s'agafa el 1-epsilon% de variabilitat
       
     if (rk!=0) eff.rank <- min(rk,Fw.r)   
     else eff.rank<-Fw.r
   }      
   
   # used rel.gvar
   used_rel.gvar<-Fw.dnou[eff.rank]
   
   # takes the effective rank dimensions of the Fw desocmpocition
   Fw.uk<-Fw.u[,1:eff.rank]  
   Fw.vk<-Fw.v[,1:eff.rank]
   Fw.dk<-Fw.d[1:eff.rank]
   
   # inverse of the singular values to prepare the Moore-Penrose pseudo-inverse of Fw.
   dFw.dk<-diag(1/Fw.dk,nrow=eff.rank,ncol=eff.rank)   

   # the Moore-Penrose pseudo-inverse of Fw and the hat matrix projector.
   Fwplus <- Fw.vk%*%dFw.dk%*%t(Fw.uk)  
   Hw<- G %*% Dsqw %*% Fwplus %*% Dsqw 
   
   # return the hat matrix, the Moore-Penrose pseudo-inverse of Fw and the effective rank.
   return(list(Hw=Hw,Fwplus=Fwplus,eff.rank=eff.rank,used_rel.gvar=used_rel.gvar))  
 }
 