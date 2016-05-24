ProjectionQuadMethod <- function(InputDim,pmaxi,Designs,Output){
  # number of PCE termrs
  M <- getM(InputDim,pmaxi) 
  # generate PCE coefficients
  PCE = generatePCEcoeff(M,Designs$Z,Output,Designs$PolyNodes,Designs$NDWi) 
  return(PCE)
}