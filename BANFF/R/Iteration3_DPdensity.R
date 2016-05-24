
Iteration3_DPdensity<-function(iter,wholeindex,dpdensitycluster,net,pirhopair,choice,rstat,v,show.steps,trace){
  z<-wholeindex
  total<-matrix(rep(0,length(wholeindex)*iter),ncol=length(wholeindex),nrow=iter)
  pro1<-0
  pro0<-0
  
  mu0=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==0)])))
  mu1=sapply(1:v, function(kk) return(mean(rstat[which(dpdensitycluster[kk,]==1)])))
  var0=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==0)])))
  var1=sapply(1:v, function(kk) return(var(rstat[which(dpdensitycluster[kk,]==1)])))
  for(jj in 1: iter){
    
    if(jj%%show.steps==0){
      cat("iter: ",jj,"\n")
      flush.console()
    }
    for(num1 in 1:length(wholeindex)){
      ztemp=c()
      
      idx = which(net[num1,]==1)
      pro0=sum((z[idx]==0))
      pro1=sum((z[idx]==1))
      
      for (kk in 1:v){
        
        if (is.na(var0[kk])){var0[kk]=var1[kk]
        }else if (is.na(var1[kk])){var1[kk]=var0[kk]}
        
        log0<-log(pirhopair$pi0[choice])+2*pirhopair$rho0[choice]*pro0-(rstat[num1]-mu0[kk])^2/(2*var0[kk])-log(sqrt(2*pi*var0[kk]))
        log1<-log(1-pirhopair$pi0[choice])+2*pirhopair$rho1[choice]*pro1-(rstat[num1]-mu1[kk])^2/(2*var1[kk])-log(sqrt(2*pi*var1[kk]))
        
        p0<-1/(1+exp(log1-log0))
        p1<-1-p0
        if (is.na(p0) | is.na(p1)) {ztemp=ztemp
        }else{
          ztemp<-c(ztemp,sample(c(0,1),1,prob=c(p0,p1)))
        }}
      if (is.null(ztemp)) {ztemp=sample(c(1,0),1,prob=c(0.5,0.5))###in case the mu is null
      }else{
        z[num1]<-sample(ztemp,1,prob=c(rep(1/length(ztemp),length(ztemp))))}
      
      
      total[jj,num1]<-z[num1]
      
      #ratio[jj,num1]<-total[jj,num1]/jj
      
    }
  }
  
  return(total)
}
