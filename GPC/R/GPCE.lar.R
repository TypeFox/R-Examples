GPCE.lar <- function(Model=NULL,
                     PCSpace="Uniform",
                     InputDim=3,
                     InputDistrib=c(),
                     ParamDistrib=NULL,
                     Q2tgt=1-10^(-6),
                     Eps=10^(-2)*(1-Q2tgt),
                     EpsForw=Eps,
                     EpsBack=Eps,
                     EnrichStep=50,
                     jmax=InputDim,
                     pmaxi=6,
                     DesignLength=8,
                     SeedSob=sample(1:1000,1)){
  
  #ishigami.fun(matrix(c(1,1,1),nrow=1))
  #Copule
  #CorreMatrix
  #DesignInput
  #OutputRe
  q=1
  #sensitivity
  ### initialization Output
  ResultObject =list()
  
  
  ### Check Distribution
  OkDistrib=c("Uniform","Gaussian","Gamma","Beta")
  if (PCSpace!="FromData"){
    if(prod(InputDistrib%in%OkDistrib)==0) {
      print("Error1: the given physic distributions cannot be traited")
      return()
    }
  }
  
  ##########
  ### Copule Formalism: To be implemented 
  ##########
  
  ### Get Design
  Args=list(Model=Model,InputDim=InputDim,InputDistrib=InputDistrib,ParamDistrib=ParamDistrib,pmaxi=pmaxi,jmax=jmax,PCSpace=PCSpace,Q2tgt=Q2tgt,SeedSob=SeedSob,EpsForw=EpsForw,EpsBack=EpsBack,EnrichStep=EnrichStep)
  ResultObject$Designs <- c(GetDesignReg(DesignLength,InputDim,SeedSob,InputDistrib,ParamDistrib,PCSpace),DesignLength=DesignLength)     
  
  x= c(ResultObject,list(Design2Eval=ResultObject$Designs$Physic,Args=Args),list(p=0,SaveQmax=c(-100,-100,-100),SaveAmax=list(c(0),c(0),c(0))))
  class(x) <- "GPCE.lar"
  
  ### Case Model=NULL return Input design
  if(is.null(Model)) { 
    return(x)
  }
  
  ### Get Output when Model is given
  Output <- Model(ResultObject$Designs$Physic)
  
  ### 
  ResultObject <- tell.GPCE.lar(x,Output)
  #Model,Args,Out$Designs,Output)
  
  return(ResultObject)
  
  
}
