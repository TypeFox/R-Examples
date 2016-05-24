RegressionLar <- function(A,DesignPCE,InputDistrib,pmaxi,Output,PCSpace,p){
  #print("Entering Regression")
  restartAnalysis<-FALSE
  DM=DataMatrix(A,DesignPCE,InputDistrib,pmaxi,PCSpace)
  IM=t(DM)%*%DM  
  #if(rcond(IM)<10^(-4)){
  #   if(kappa(IM)>10^(5)){    
  #     R2=-10
  #     Q2=-10
  #     restartAnalysis<-TRUE
  #     return(list(R2=R2,Q2=Q2,rAnaly=restartAnalysis))
  #   } 
  #   else {
  #print(IM)
  InvIM=solve(IM)
  a=InvIM%*%t(DM)%*%t(t(Output))
  MetaOut=t(a)%*%t(DM) #ligne
  R2=EE(MetaOut,Output)
  #R2adj=EEadj(MetaOut)
  #     print(DM)
  #     print(IM)
  #     print(diag(DM%*%InvIM%*%t(DM)))
  h=diag(DM%*%InvIM%*%t(DM))
  if (1 %in% h) {print('h problem'); Q2=-100; restartAnalysis=TRUE}
  else {Q2=LOOLar(MetaOut,diag(DM%*%InvIM%*%t(DM)),Output,InvIM,length(InputDistrib),p)}
  return(list(R2=R2,Q2=Q2,rAnaly=restartAnalysis))
  #   }  
}