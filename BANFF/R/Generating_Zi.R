#####Generating Zi by different setting of pi and rho
Generating_Zi=function(net,pirhopair,wholeindex,n=30)
{
  z=wholeindex
  znew=c()
  
  for (j in 1:length(pirhopair$pi0))
  {
    for(i in 1: n)
      
    {
      pro1<-0
      pro0<-0
      for(num1 in 1:length(wholeindex)){
        idx = which(net[num1,]==1)
        pro0=sum(z[idx]==0)
        pro1=sum(z[idx]==1)
        
        
        log0<-log(pirhopair$pi0[j])+(2*pirhopair$rho0[j]*pro0)
        
        log1<-log(1-pirhopair$pi0[j])+(2*pirhopair$rho1[j]*pro1)
        
        
        
        p0<-1/(1+exp(log1-log0))
        p1<-1-p0
        
        z[num1]<-sample(c(0,1),1,prob=c(p0,p1))
      }
    }
    znew=rbind(znew,z)
  }
  return(znew)
}  

