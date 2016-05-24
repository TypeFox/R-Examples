mc.spGini<-function(Nsim=99,Bandwidth,x,Coord.X,Coord.Y,WType='Binary')
{
  SIM.m<-matrix(data=NA, nrow=Nsim, ncol=5)
  
  Obs<-length(x)
  
  Coords<-cbind(Coord.X,Coord.Y)
  
  spGini.Obs<-spGini(Coords,Bandwidth,x,WType)
  
  SG<-spGini.Obs[[3]]/spGini.Obs[[1]]
  
  for(i in 1:Nsim){ 
   
    spGini.SIM<-spGini(Coords[sample(nrow(Coords)),],Bandwidth,x,WType)
        
    SIM.m[i,1]<-i
    SIM.m[i,2]<-spGini.SIM[[2]]
    SIM.m[i,3]<-spGini.SIM[[3]]
    SIM.m[i,4]<-spGini.SIM[[3]] / spGini.SIM[[1]]
    if (SIM.m[i,4]>=SG){
      SIM.m[i,5]<-1}
    else
    {SIM.m[i,5]<-0}
}
  C.test<-sum(SIM.m[,5])
  
  if ((Nsim-C.test)< C.test) {C.test=Nsim-C.test}
  
  pseudo.p<-(1+C.test)/(Nsim+1)
  
  SIMs<-data.frame(SIM.ID=SIM.m[,1], SIM.gwGini= SIM.m[,2],SIM.nsGini=SIM.m[,3], SIM.SG=SIM.m[,4], SIM.Extr= SIM.m[,5])
  
  return(list(SIM=SIMs,spGini.Observed=spGini.Obs,pseudo.p=pseudo.p))
}