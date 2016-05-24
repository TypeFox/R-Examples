#####################################################################################################################
#################### Find optimal delta.star for given lambda #######################################################
################## Reference: SLedoit and Wolf (2003) and Sch\"afer and Strimmer (2005) #############################


delta.star<-function(data,thres.range)
  
{
    
  n<-dim(data)[1]
  
  p<-dim(data)[2]
  
  cov.hat<-cov(data)
  
  mean<-colMeans(data)
  
  stanX<-scale(data)  # standardized data matrix
  
  cov.soft<-softt(cov(stanX),thres.range)
  
  w<-array(0,c(p,p,n))
  
  w.mean<-Var.S<-matrix(0,p,p)
  
  l<-length(thres.range)
  
  delta.star<-Nu<-De<-rep(0,l)
  
  cov.novel<-I<-array(cov.hat,c(p,p,l))
  
  cor.novel<-array(cov2cor(cov.hat),c(p,p,l))  
  
  for (i in 1:p)
  
  {
   
    for (j in i:p)
      
    {
      
      w[i,j,]<- w[j,i,]<-stanX[,i]*stanX[,j]
      
      w.mean[i,j]<-w.mean[j,i]<-mean(w[i,j,])
   
    }
    
  }
  
  S<-n/(n-1)*w.mean  # as the same as S<-cov(X)
  
  for (i in 1:p)
    
  {
    
    for (j in i:p)
      
    {
      
      Var.S[i,j]<-Var.S[j,i]<-n/(n-1)^3*sum((w[i,j,]-w.mean[i,j])^2)
      
    }
    
  }
  
  VarS.hat<-Var.S-diag(diag(Var.S))
  
  for (i in 2:l)
 
  {
    
    I[,,i]<-abs(S)<thres.range[i]
    
    Nu[i]<-sum(VarS.hat*I[,,i])
    
    De[i]<-sum((S-cov.soft[,,i])^2)
    
    delta.star[i]<-Nu[i]/De[i]
        
    cor.novel[,,i]=(1-delta.star[i])*S+delta.star[i]*(cov.soft[,,i])
    
    cov.novel[,,i]=sqrt(diag(diag(cov.hat)))%*%cor.novel[,,i]%*%sqrt(diag(diag(cov.hat)))
    
  }
  
  th.del.star.cov<-matrix(c(thres.range,delta.star),l,2)
  	
  colnames(th.del.star.cov)<-c("thres","delta.star")

  return(list("covariance.novelist.candidates"=cov.novel,"correlation.novelist.candidates"=cor.novel,"delta.star"=th.del.star.cov))
  
}

