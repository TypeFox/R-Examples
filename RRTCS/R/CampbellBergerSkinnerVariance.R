CampbellBergerSkinnerVariance=function(N,r,pi,pij,type=c("total","mean")){
  if(type=="total"){
   return(VE.Jk.CBS.SYG.Total.Hajek(r,pi,pij,N))
  }
  if(type=="mean"){
    return(VE.Jk.CBS.SYG.Mean.Hajek(r,pi,pij))
  }
}