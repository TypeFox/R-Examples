PostProcessProj <- function(Out,InputDim,pmaxi){

  ### Number of full PCE
  ### Mean, Variance and stand. Dev
  M <- getM(InputDim,pmaxi)
  
  PCEMean = Out$PCEcoeff[1]
  PCEVar = sum(Out$PCEcoeff[2:M]^2*Out$PhiIJ[2:M])
  PCESD = sqrt(PCEVar)
  
  #Skew=mean((MOut-PCEMean)^3)/PCESD^3
  #Kurto=mean((MOut-PCEMean)^4)/PCESD^4
    
  # Sobol' sensitivity analysis
  Index = indexCardinal(InputDim,pmaxi)
  Sobol = getSobol(InputDim,Index,Out$PCEcoeff,Out$PhiIJ)
  
  ### Output distribution
  #H.pi <- dpik(MOut)
  #k.de <- bkde(MOut,bandwidth=H.pi)
  return(list(Moments=list(PCEMean=PCEMean,PCEVar=PCEVar,PCESD=PCESD),Sensitivity=Sobol))
}