novelist.assign.inv<- function(m, th, delta,thf = softt) 
  
{
  
  # NOVELIST with change of thresholding and shrinkage intensity as well
  # th is a vector of thresholding level
  # delta is a vector of shrinkage intensity 
  
  diag.m<-diag(diag(m))
  
  m<-cov2cor(m)
  
  l1<-length(delta)
  
  l <- length(th)
  
  dm <- dim(m)

  p<-dm[1]
  
  novel.cor.all <- novel.cov.all <-novel.cov.inv.all<-novel.cor.inv.all<- array(matrix(NA,p, p), c(p, p, l,l1))
  
  for (i in 1:l) 
  
  {
    
    for (j in 1:l1)
    
    {

     
      novel.cor.all[,,i,j]<-(1-delta[j]) * m+ delta[j] * thf(m, th)[,,i]
     
      novel.cov.all[,,i,j]<-sqrt(diag.m)%*%novel.cor.all[,,i,j]%*%sqrt(diag.m)

      e<-eigen(novel.cov.all[,,i,j])$values
      
      if(e[p]>0&((e[1]/e[p])<10^6)) 

      {
      
      novel.cov.inv.all[,,i,j]<-solve(novel.cov.all[,,i,j])

      }

      e<-eigen(novel.cor.all[,,i,j])$values
      
      if(e[p]>0&(e[1]/e[p])<10^6) 

      {
      
      novel.cor.inv.all[,,i,j]<-solve(novel.cor.all[,,i,j])

      }

 
    }
    
  }
  
  return(list("inv.cor.novel"= novel.cor.inv.all,"inv.cov.novel"= novel.cov.inv.all))
  
}