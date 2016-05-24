RegressionSparseMethod <- function(InputDim,Designs,InputDistrib,ParamDistrib,pmaxi,jmax,Output,PCSpace,
                                   Q2tgt,SeedSob,EpsForw,EpsBack,EnrichStep){
  #print("Entering RegressionSparseMethod")
  #print(Designs$DesignLength)
  restartAnalysis<-FALSE    
  p = 0
  Ap = matrix(rep(0,InputDim),1,InputDim)
  A = Ap
  Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
  R2Q2_0 <- c(Reg$R2,Reg$Q2)
  #print(R2Q2_0)
  #print(p)
  while( R2Q2_0[2]<=Q2tgt & p<pmaxi){      
    p = p+1
    #print(c("poly degree:",p))
    j = 1
    App=c()
    while( R2Q2_0[2]<=Q2tgt & j<=min(p,jmax)){
      #print(c("interaction:", j))
      J_jp=c()
      TC_jp = HIS(p,InputDim,1)#t(indexCardinal(InputDim,p))
      C_jp = matrix(TC_jp[which( rowSums(TC_jp>0)==j & rowSums(TC_jp)==p ),],ncol=InputDim)
      k=0
      for(i in 1:nrow(C_jp)){
        A = rbind(Ap,C_jp[i,])
        #print(c("p:",p,"j:",j,"C_jp:",C_jp[i,]))
        Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
        R2Q2 <- c(Reg$R2,Reg$Q2)
        restartAnalysis <- Reg$rAnaly
        if(restartAnalysis==TRUE){
          #print("Enrich2");
          Designs$DesignLength <- Designs$DesignLength + EnrichStep
          ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)
          return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength)))
        }
        DeltaR2 = R2Q2[1]-R2Q2_0[1]
        if(DeltaR2>=EpsForw){ J_jp = rbind(J_jp,c(C_jp[i,],DeltaR2)); k=k+1 }          
      }
      if(k>0 & restartAnalysis==FALSE){
        J_jpStar=matrix(J_jp[sort(J_jp[,InputDim+1],index.return=TRUE,decreasing=TRUE)$ix,(1:InputDim)],ncol=InputDim)
        R_jp=c()
        for(i in 1:nrow(J_jpStar)){
          #print(c("p:",p,"j:",j,"jpStar:",J_jpStar[i,]))
          DM=DataMatrix(rbind(Ap,J_jpStar[i,]),Designs$PCE,InputDistrib,pmaxi,PCSpace)
          IM=t(DM)%*%DM
          if(rcond(IM)<10^(-4)){
            #print("Enrich1")
            Designs$DesignLength <- Designs$DesignLength + EnrichStep
            ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)
            return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength)))
          }
          else { 
            R_jp = rbind(R_jp,J_jpStar[i,]) 
          }
        }
        App = rbind(Ap,R_jp)
        Ap=App
      }
      j=j+1        
    }
    if(restartAnalysis==FALSE){
      Reg <- Regression(Ap,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
      R2Q2_0 <- c(Reg$R2,Reg$Q2)
      restartAnalysis <- Reg$rAnaly
      IndSup=c()
      for(i in 1:nrow(Ap)){
        #print(c("p:",p,"j:",j,"Ap minus",Ap[i,]))
        A=matrix(Ap[-i,],ncol=InputDim)
        Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
        R2Q2 <- c(Reg$R2,Reg$Q2)
        restartAnalysis <- Reg$rAnaly
        if(restartAnalysis==TRUE){
          #print("Enrich3");
          Designs$DesignLength <- Designs$DesignLength + EnrichStep
          ED<-GetDesignReg(Designs$DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace)
          return(list(Designs=list(PCE=ED$PCE,Physic=ED$Physic,Sobol=ED$Sobol,DesignLength=Designs$DesignLength)))
        }
        DeltaR2 = R2Q2_0[1]-R2Q2[1]
        if(DeltaR2<=EpsBack){ IndSup=cbind(IndSup,i) }
      }
      if(length(IndSup)!=0 & restartAnalysis==FALSE){
        Ap<-Ap[-IndSup,]
        Reg <- Regression(A,Designs$PCE,InputDistrib,pmaxi,Output,PCSpace)
        R2Q2_0 <- c(Reg$R2,Reg$Q2)
        restartAnalysis <- Reg$rAnaly
      }
    }      
  }
  #print("Leaving RegressionSparse...")
  return(list(TruncSet=Ap,R2=R2Q2_0[1],Q2=R2Q2_0[2],Designs=list(PCE=Designs$PCE,Physic=Designs$Physic,Sobol=Designs$Sobol,DesignLength=Designs$DesignLength,p=p)))    
}  

