EscobarBergerVariance=function(N,r,pi,pij,type=c("total","mean")){
  if(type=="total"){
    return(VE.EB.SYG.Total.Hajek(r,pi,pij,N,rep(1,length(pi))))
  }
  if(type=="mean"){
    return(VE.EB.SYG.Mean.Hajek(r,pi,pij,rep(1,length(pi))))
  }
}