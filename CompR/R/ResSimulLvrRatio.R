ResSimulLvrRatio<-function(Data,ResH0,Constraint=0,nsimul,level=0.05,eps=0.0001,eps1=0.0001)
  
{
  ResLvrRatio<-new(Class="LvrRatio")
  TclaH0<-length(ResH0@Pi)
  ResLvrRatio@Simu<-matrix(0,nsimul,4)
  colnames(ResLvrRatio@Simu)<-c("TclaH0","LvrH0","LvrH1","LvrH1-LvrH0")
  ResLvrRatio@Test<-matrix(0,1,2)
  colnames(ResLvrRatio@Test)<-c("Level","Quantile")
  for (simul in 1:nsimul)
  {
    DataSimul<-DataSimulH0(Data,ResH0)
    Resglobsimu0<-EstimBradley(DataSimul,Constraint,TclaH0,TestPi=FALSE,eps=eps,eps1=eps1)
    Resglobsimu1<-EstimBradley(DataSimul,Constraint,TclaH0+1,TestPi=FALSE,eps=eps,eps1=eps1)
    ResLvrRatio@Simu[simul,1]<-TclaH0
    ResLvrRatio@Simu[simul,2]<-Resglobsimu0@Lvr
    ResLvrRatio@Simu[simul,3]<-Resglobsimu1@Lvr
  }
  ResLvrRatio@Simu[,4]<-ResLvrRatio@Simu[,3]-ResLvrRatio@Simu[,2]
  ResLvrRatio@Test[1,1]<-level
  ResLvrRatio@Test[1,2]<-quantile(ResLvrRatio@Simu[,4],probs=1-level)
  
  return(ResLvrRatio)
}
