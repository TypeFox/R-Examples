predict.GPCE.lar <- function(object,x,Selection=1:nrow(object$TruncSet),...) {
  A=object$TruncSet[Selection,]
  Design=GetDesign4Predict(x,object$Args$InputDistrib,object$Args$ParamDistrib,object$Args$PCSpace)
  DM=DataMatrix(A,Design$PCE,object$Args$InputDistrib,object$Args$pmaxi,object$Args$PCSpace)
  MetaOut=DM%*%object$CoeffPCE[Selection,]  
  return(MetaOut)
}
