PostProcessReg <- function(A,DesignPCE,InputDim,InputDistrib,ParamDistrib,pmaxi,Output,jmax,PCSpace){
  print("Entering PostProcessReg")
  ### Number of PC
  #P <- factorial(InputDim + pmaxi)/(factorial(InputDim)*factorial(pmaxi))
  ### Final Decomposition Coefficients
  DM=DataMatrix(A,DesignPCE,InputDistrib,pmaxi,PCSpace)
  IM=t(DM)%*%DM
  InvIM=solve(IM)
  a=InvIM%*%t(DM)%*%t(t(Output))
  
  ### Mean, Variance and stand. Dev
  #MeanO=mean(Output)
  #VarO=DesignLength*var(Output)/(DesignLength-1)
  #SDO=sqrt(VarO)
  PCEMean=a[1]
  PCEVar=sum(a[-1]^2)
  PCESD=sqrt(PCEVar)
  
  ### Skewness and Kurtosis
  NDesignLength=5000
  NDesign=GetDesignReg(NDesignLength,InputDim,SeedSob=sample(1:10000, 1),InputDistrib,ParamDistrib,PCSpace)
  MOut=EvalPCE(a,A,NDesign$PCE,InputDistrib,pmaxi,PCSpace)  
  Skew=mean((MOut-PCEMean)^3)/PCESD^3
  Kurto=mean((MOut-PCEMean)^4)/PCESD^4
  ### Output distribution
  #H.pi <- dpik(MOut)
  #k.de <- bkde(MOut,bandwidth=H.pi)
  H.pi=0
  k.de=0
  ### Sensitivity indices
  S=list()
  for(i in 1:jmax){
    Comb=combn(InputDim,i)
    tempS=matrix(0,nrow(Comb)+1,ncol(Comb))
    tempS[-1,]=Comb
    for (j in 1: ncol(Comb)){    
      Ind=which( apply(as.matrix(A[,-Comb[,j]]),1,FUN="sum")==0 & apply(as.matrix(A[,Comb[,j]]),1,FUN="prod")>0)    
      if(length(Ind)!=0) tempS[1,j]=sum(a[Ind]^2)/PCEVar    
    }
    S=c(S,list(tempS))
  }
  
  ### Total Sensitivity Indices
  ST=list()
  for(i in 1:jmax){
    Comb=combn(InputDim,i)
    tempS=matrix(0,nrow(Comb)+1,ncol(Comb))
    tempS[-1,]=Comb
    for (j in 1: ncol(Comb)){    
      Ind=which( apply(as.matrix(A[,Comb[,j]]),1,FUN="prod")>0)    
      if(length(Ind)!=0) tempS[1,j]=sum(a[Ind]^2)/PCEVar    
    }
    ST=c(ST,list(tempS))
  }
  return(list(CoeffPCE=a,Moments=list(PCEMean=PCEMean,PCEVar=PCEVar,PCESD=PCESD,Skew=Skew,Kurto=Kurto),Sensitivity=list(S=S,ST=ST),OutputDistib=list(bandw=H.pi,Distr=k.de)))
}