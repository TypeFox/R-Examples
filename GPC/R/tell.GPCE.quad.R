tell.GPCE.quad <- function(x,Output,...){ #ResultObjectDesign2EvalArgs,Output){  
  ResultObjectDesign2EvalArgs=x

  Args    = ResultObjectDesign2EvalArgs$Args
  Design  = ResultObjectDesign2EvalArgs$Design
  Model   = ResultObjectDesign2EvalArgs$Args$Model
  
  ### Expansion Coeff estimation
  ### Projection method with Quadrature scheme
  M <- getM(Args$InputDim,Args$p)
  ResultObject <- generatePCEcoeff(M,Design$QuadSize,Output,Design$PolyNodes,Design$QuadWeights)
  
  if(names(ResultObject)[1]=="PCEcoeff") { # return results
    OutProcess <- PostProcessProj(ResultObject,Args$InputDim,Args$p)     
    x0 <-c(ResultObject,OutProcess,list(Args=Args),list(Output=Output),list(Design=Design))
    class(x0) <- "GPCE.quad"
    return(x0)
  }
  if(is.null(Model)) { 
    x1=c(ResultObject,list(Design2Eval=ResultObject$Design$PhysicalNodes[(ResultObject$Design$DesignLength-Args$EnrichStep+1):ResultObject$Designs$DesignLength,],Output=Output,Args=Args)) 
    class(x1) <- "GPCE.quad"
    return(x1)
  }
}

