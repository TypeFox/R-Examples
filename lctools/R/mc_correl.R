mc.lcorrel <- function(Nsim=99,bwSIM,CorVars,Coord.X,Coord.Y) 
{

  lcorSIM<-matrix(data=NA, nrow=Nsim, ncol=11)
  
  Obs<-nrow(CorVars)
  
  Coords<-cbind(Coord.X,Coord.Y)

  lcc0<-lcorrel(CorVars,bwSIM, Coords)
  Obs.var.lpcc<-var(lcc0$LPCC[,2])
  Obs.var.gwpcc<-var(lcc0$GWPCC[,2])

  for(i in 1:Nsim){ 
    
    lccSIM<-lcorrel(CorVars,bwSIM,Coords[sample(nrow(Coords)),])
    
    lcorSIM[i,1]<-i
    lcorSIM[i,2]<-mean(lccSIM$LPCC[,2])
    lcorSIM[i,3]<-min(lccSIM$LPCC[,2])
    lcorSIM[i,4]<-max(lccSIM$LPCC[,2])
    lcorSIM[i,5]<-var(lccSIM$LPCC[,2])
    
    if (lcorSIM[i,5]>=Obs.var.lpcc){
      lcorSIM[i,6]<-1}
    else
    {lcorSIM[i,6]<-0}
    
    lcorSIM[i,7]<-mean(lccSIM$GWPCC[,2])
    lcorSIM[i,8]<-min(lccSIM$GWPCC[,2])
    lcorSIM[i,9]<-max(lccSIM$GWPCC[,2])  
    lcorSIM[i,10]<-var(lccSIM$GWPCC[,2]) 
    
    if (lcorSIM[i,10]>=Obs.var.gwpcc){
      lcorSIM[i,11]<-1}
    else
    {lcorSIM[i,11]<-0}

  }
  
  C.test.lp<-sum(lcorSIM[,6])
  C.test.gw<-sum(lcorSIM[,11])
  
  if ((Nsim-C.test.lp)< C.test.lp) {C.test.lp=Nsim-C.test.lp}
  if ((Nsim-C.test.gw)< C.test.gw) {C.test.gw=Nsim-C.test.gw}
  
  pseudo.p.lpcc = (1+C.test.lp)/(Nsim+1)
  pseudo.p.gwpcc = (1+C.test.gw)/(Nsim+1)
  
  SIMs<-data.frame(SIM.ID=lcorSIM[,1], SIM.meanLPCC=lcorSIM[,2],SIM.minLPCC=lcorSIM[,3], SIM.maxLPCC=lcorSIM[,4],SIM.varLPCC=lcorSIM[,5], SIM.extrLPCC= lcorSIM[,6], SIM.meanGWPCC=lcorSIM[,7],SIM.minGWPCC=lcorSIM[,8], SIM.maxGWPCC=lcorSIM[,9],SIM.varGWPCC=lcorSIM[,10], SIM.extrGWPCC= lcorSIM[,11])
    
  return(list(SIM=SIMs,LC.Obs=lcc0, pseudo.p.lpcc = pseudo.p.lpcc,pseudo.p.gwpcc =pseudo.p.gwpcc ))
}

  

