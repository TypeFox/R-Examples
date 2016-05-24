ce.simnormal.Init.Mean.AIC <-
function(N, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init){
    
    v <- var(data[, 1])
    
    #   if (N==0){
    #     seql<-c(1,L)
    #     mBic.full<-mBIC(seql,data,0,L,h)
    #     
    #     return(list(locis=c(1,L+1),mBIC=mBic.full))
    #     rm(mBic.full,seql)
    #     
    #   } else {
    
    ########################Parameter initialization######################################################  
    #new_para<-rbind(rep(L0+(L-L0)/2,N),rep(sqrt(L-L0)^2/12,N))    
    new_para <- rbind(init.locs, rep(var.init, N))    
    
    #  n_par_m <- array(init.locs, dim=c(1,N))   
    #  n_par_sd<-array(std,dim=c(1,N))   
    #  new_para<-rbind(n_par_m,n_par_sd)   
    ######################################################################################################  
#     llVal <- c()
#     aic <- c()
    k<-0
    repeat
    {
      k<-k+1
      ch<-array(0,dim=c(M,N+2))       
      ch[,1]<-c(1)                    
      ch[,N+2]<-c(L+1)    
      ch[,(2:(N+1))]<-apply(new_para,2,normrand,L0,L,M)       
      ch<-t(apply(ch,1,sort))    
      
#       LL.full <- apply(ch, 1, loglikMeanNormal, data, h)
#       AIC.val <- -2*LL.full + 2*(N + 2)
      LL.full <- apply(ch, 1, llhood.MeanNormal, data, v, h)
      AIC.val <- apply(as.data.frame(LL.full), 1, AIC.MeanNormal, N, L) 
      
      #mod_bic<-apply(ch,1,mBIC,data,N,L,h)
      #ch<-cbind(ch,mod_bic)                         
      ch <- cbind(ch, LL.full, AIC.val)
      ch <- ch[order(ch[, (N + 4)], decreasing = FALSE), ]  
      
      melitesmpl<-ch[1:Melite,]                     
#       llVal[k] <- melitesmpl[1, (N + 3)] 
#       aic[k] <- melitesmpl[1, (N + 4)] 
      
      
      new_par_n<-array(0,dim=c(2,N))
      new_par_n[1,]<-apply(as.matrix(melitesmpl[,(2:(N+1))]),2,mean)
      new_par_n[2,]<-apply(as.matrix(melitesmpl[,(2:(N+1))]),2,sd)   
      
      new_para[1,] <- a*new_par_n[1,] + (1-a)*new_para[1,]
      new_para[2,] <- b*new_par_n[2,] + (1-b)*new_para[2,]
      
      mad<-apply(as.matrix(melitesmpl[,(2:(N+1))]),2,mad)
      
      if(max(mad)<=eps){break}
    }
#    return(list(loci=ch[1,(1:(N+2))], AIC.Val = aic[k], LogLike = llVal[k]))
    return(list(loci=ch[1,(1:(N+2))], AIC.Val = melitesmpl[1, (N + 4)][[1]], LogLike = melitesmpl[1, (N + 3)][[1]]))
    #  }
  }
