EvalPCE <- function(a,A,Design,InputDistrib,pmaxi,PCSpace){
  print("Entering EvalPCE")
  DM=DataMatrix(A,Design,InputDistrib,pmaxi,PCSpace)
  MetaOut=DM%*%a #ligne
  return(MetaOut)
}