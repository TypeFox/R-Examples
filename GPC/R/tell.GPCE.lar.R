tell.GPCE.lar <- function(x,Output,...){ #ResultObjectDesign2EvalArgs,Output){
  ResultObjectDesign2EvalArgs=x
  
  Args = ResultObjectDesign2EvalArgs$Args
  Designs = ResultObjectDesign2EvalArgs$Designs
  Model = Args$Model
  
  ### Expansion Coeff estimation
  ### Regression lar expansion method
  ResultObject <- RegressionlarMethod(Args$InputDim,Designs,Args$InputDistrib,Args$ParamDistrib,
                                      Args$pmaxi,Args$jmax,Output,Args$PCSpace,Args$Q2tgt,Args$SeedSob,Args$EpsForw,
                                      Args$EpsBack,Args$EnrichStep,ResultObjectDesign2EvalArgs$p,
                                      ResultObjectDesign2EvalArgs$SaveQmax,ResultObjectDesign2EvalArgs$SaveAmax)
  if(names(ResultObject)[1]=="TruncSet") {
    OutProcess <- PostProcessReg(ResultObject$TruncSet,ResultObject$Designs$PCE,Args$InputDim,Args$InputDistrib,Args$ParamDistrib,Args$pmaxi,Output,Args$jmax,Args$PCSpace)     
    x0 <-c(ResultObject,OutProcess,list(Args=Args),list(Output=Output))
    class(x0) <- "GPCE.lar"
    return(x0)
  }
  if(is.null(Model)) { 
    x1=c(ResultObject,
         list(Design2Eval=ResultObject$Designs$Physic[(ResultObject$Designs$DesignLength-Args$EnrichStep+1):ResultObject$Designs$DesignLength,],
              Output=Output,Args=Args)) 
    #c(Out,list(Design2Eval=Out$Designs$Physic,Args=Args))
    class(x1) <- "GPCE.lar"
    return(x1)
  }
  else { 
    x= c(ResultObject,list(Design2Eval=ResultObject$Designs$Physic,Args=Args))
    class(x) <- "GPCE.lar"
    Output <- c(Output,Model(ResultObject$Designs$Physic[(ResultObject$Designs$DesignLength-Args$EnrichStep+1):ResultObject$Designs$DesignLength,]))
    tell.GPCE.lar(x,Output)
  }
}