novelist.assign<- function(m, th, delta,thf = softt) 
  
{
  
  # NOVELIST with change of thresholding and shrinkage intensity as well
  # th is a vector of thresholding level
  # delta is a vector of shrinkage intensity 
  
  diag.m<-diag(diag(m))
  
  m<-cov2cor(m)
  
  l1<-length(delta)
  
  l <- length(th)
  
  dm <- dim(m)
  
  novel.cor.all <- novel.cov.all <- array(m, c(dm[1], dm[2], l,l1))
  
  for (i in 1:l) 
  
  {
    
    for (j in 1:l1)
    
    {


      novel.cor.all[,,i,j]<-(1-delta[j]) * m+ delta[j] * thf(m, th)[,,i]
     
      novel.cov.all[,,i,j]<-sqrt(diag.m)%*%novel.cor.all[,,i,j]%*%sqrt(diag.m)


 
    }
    
  }
  
  return(list("cor.novel"= novel.cor.all,"cov.novel"= novel.cov.all))
  
}