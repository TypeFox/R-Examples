GetDesignReg <- function(DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace){
  #print("Entering GetDesignReg")
  ### Design Sobol [0,1]
  DesignSobol <- randtoolbox::sobol(n=DesignLength, dim = InputDim, init = TRUE, scrambling = 1, seed = SeedSob, normal = FALSE)
  ### Init DesignPhysic 
  DesignPhysic <- DesignSobol 
  
  ### Get DesignPhysic
  for (i in 1:InputDim){
    if (InputDistrib[i]=="Gaussian"){
      DesignPhysic[,i]=qnorm(DesignSobol[,i],0,1)
    }
    if (InputDistrib[i]=="Gamma"){
      DesignPhysic[,i]=qgamma(DesignSobol[,i],shape=ParamDistrib$shapeG,scale=ParamDistrib$scaleG)
    }
    if (InputDistrib[i]=="Beta"){
      DesignPhysic[,i]=qbeta(DesignSobol[,i],shape1=ParamDistrib$shape1B,shape2=ParamDistrib$shapeB)
    } 
  }
  ### Get Design for PCE expansion
  if(PCSpace=="Uniform"){ DesignPCE <- 2*DesignSobol-1 }
  if(PCSpace=="Gaussian"){ DesignPCE <- qnorm(DesignSobol,0,1) }
  if(PCSpace=="Physic"){ DesignPCE <- DesignPhysic }
  return(list(Physic=DesignPhysic,PCE=DesignPCE,Sobol=DesignSobol))
}
