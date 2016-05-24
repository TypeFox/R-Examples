ce.simnormal.Init.mBIC <-
function(N, init.locs, data, h, L0, L, M, Melite, eps, a, b, var.init){
  
  if (N==0){
    seql<-c(1,L)
    mBic.full<-mBIC(seql,data,0,L,h)
    
    return(list(locis=c(1,L+1),mBIC=mBic.full))
    rm(mBic.full,seql)
    
  } else {
    
  ########################Parameter initialization######################################################  
    new_para <- rbind(init.locs, rep(var.init, N))    
  ######################################################################################################  
#  modbic<-c()
  
  k<-0
  repeat
  {
    k<-k+1
    ch<-array(0,dim=c(M,N+2))       
    ch[,1]<-c(1)                    
    ch[,N+2]<-c(L+1)    
    ch[,(2:(N+1))]<-apply(new_para,2,normrand,L0,L,M)       
    ch<-t(apply(ch,1,sort))                       
    mod_bic<-apply(ch,1,mBIC,data,N,L,h)
    ch<-cbind(ch,mod_bic)                         
    ch<-ch[order(ch[,(N+3)],decreasing=TRUE),]  
    melitesmpl<-ch[1:Melite,]                     
#    modbic[k]<-melitesmpl[1,(N+3)] 
    
    new_par_n<-array(0,dim=c(2,N))
    new_par_n[1,]<-apply(as.matrix(melitesmpl[,(2:(N+1))]),2,mean)
    new_par_n[2,]<-apply(as.matrix(melitesmpl[,(2:(N+1))]),2,sd)   
    
    new_para[1,] <- a*new_par_n[1,] + (1-a)*new_para[1,]
    new_para[2,] <- b*new_par_n[2,] + (1-b)*new_para[2,]
    
    mad<-apply(as.matrix(melitesmpl[,(2:(N+1))]),2,mad)
    
    if(max(mad)<=eps){break}
  }
#  return(list(loci=ch[1,(1:(N+2))], mBIC=modbic[k]))
  return(list(loci=ch[1,(1:(N+2))], mBIC=melitesmpl[1,(N+3)][[1]]))
  }
}
