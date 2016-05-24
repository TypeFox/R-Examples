GetDesign4Predict <- function(Design,InputDistrib,ParamDistrib,PCSpace){
  DesignSobol <- Design
  DesignPCE <- DesignSobol 
  for (i in 1: ncol(Design)){
    if (InputDistrib[i]=="Gaussian"){
      DesignSobol[,i]=pnorm(Design[,i],0,1)
    }
    if (InputDistrib[i]=="Gamma"){
      DesignSobol[,i]=pgamma(Design[,i],shape=ParamDistrib$shapeG,scale=ParamDistrib$scaleG)
    }
    if (InputDistrib[i]=="Beta"){
      DesignSobol[,i]=pbeta(Design[,i],shape1=ParamDistrib$shape1B,shape2=ParamDistrib$shapeB)
    } 
  }
  ### Get Design for PCE expansion
  if(PCSpace=="Uniform"){ DesignPCE <- 2*DesignSobol-1 }
  if(PCSpace=="Gaussian"){ DesignPCE <- qnorm(DesignSobol,0,1) }
  if(PCSpace=="Physic"){ DesignPCE <- Design }
  return(list(Physic=Design,PCE=DesignPCE,Sobol=DesignSobol))
}
