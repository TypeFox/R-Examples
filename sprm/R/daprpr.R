daprpr <-
function(Data,center.type,scale.type){
  
  # DAPRPR Data preprocessing
  # Classical and robust scaling 
  # Inputs: Data, type of centering as R function name (e.g. "mean", "l1median")
  #         type of scaling as R function name (e.g. "sd", "qn","Sn","scaleTau2")
  # Robust estimators: Please use library(robustbase)
  # Outputs: Scaled data, attributes: Data Center, Data Scale
  # Given center.type="mean" and scale.type="sd", dpp(...) is 
  # equivalent to scale(Data, center=T, scale=T). 
  # Written by Sven Serneels, BASF, GVM/S, June 2013.
  
#  require(robustbase)
#  require(pcaPP)
  if(any(is.na(Data))==TRUE){stop("Routine cannot process NAs")}
  if(missing(center.type)){center.type="median"}
  if(missing(scale.type)){scale.type="qn"}
  if(center.type=="l1median"){
    Data.center <- l1median(Data)
  } else {
    Data.center <- apply(Data,2,center.type)
  }
  Data.scaled <- (Data - matrix(Data.center,nrow=dim(Data)[1],ncol=dim(Data)[2],byrow=TRUE))
  if(!(scale.type=="no")){
    Data.scale <- apply(Data,2,scale.type)
    if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
      Data.scale <- apply(Data,2,sd)
      if (any(1/Data.scale>1e19) | any(is.nan(Data.scale))){
      Data.scale <- rep(1,ncol(Data))
      warning("Routine used scale.type='no' to avoide division by zero or infinity.")
      } else {
        warning("Routine used scale.type='sd' to avoide division by zero or infinity.")
      }
    }
    Data.scaled <- Data.scaled/matrix(Data.scale,nrow=dim(Data)[1],ncol=dim(Data)[2],byrow=TRUE)
  } else {
    Data.scale <- rep(1,ncol(Data))
  }
  attr(Data.scaled,"Center") <- Data.center
  attr(Data.scaled,"Scale") <- Data.scale
  attr(Data.scaled,"Type") <- c(center.type,scale.type)
  return(Data.scaled)
}
