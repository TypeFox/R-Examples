DataSimulH0<-function(Data,ResH0)
{
  nsujets<-length(Data@Cons)
  ncrit<-length(Data@Crit)
  nprod<-length(Data@Prod)
  Tcla<-length(ResH0@Pi)
  DataSimulH0<-Data
  LambdaH0<-ResH0@Lambda
  Pi<-ResH0@Pi
  
  
  
  for (h in 1:nsujets)
  {
    probnew<-LambdaH0[1]
    alea<-runif(1,0,1)
    t<-1
    while (alea>probnew)
    {
      t<-t+1
      probnew<-probnew+LambdaH0[t]
    }
    Pit<-Pi[[t]]
    Matreftcrit<-contrindiv(Pit)
    for (crit in 1:ncrit)
    {
      DataSimulH0@Paircomp[[crit]][[h]]<-0*DataSimulH0@Paircomp[[crit]][[h]]
      Refhcrit<-Data@Paircomp[[crit]][[h]]+t(Data@Paircomp[[crit]][[h]])
      while (any(Refhcrit>0))
      {
        Msimulhcrit<-matrix(runif(nprod*nprod),nrow=nprod,ncol=nprod)
        Msimulhcrit<-(matrix(1,ncol=nprod,nrow=nprod))*(Msimulhcrit<Matreftcrit[[crit]])
        a<-Msimulhcrit*lower.tri(Msimulhcrit)
        b<-1-t(a)
        b<-b*upper.tri(b)
        Msimulhcrit<-a+b
        Pairviewhcrit=Refhcrit*(Refhcrit>0)
        DataSimulH0@Paircomp[[crit]][[h]]<-DataSimulH0@Paircomp[[crit]][[h]]+Msimulhcrit*Pairviewhcrit
        Refhcrit<-Refhcrit-1
      }
    }
  }
  return(DataSimulH0)
}