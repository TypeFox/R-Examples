#'
#' This function transforms the two matrices CN and fracB in one matrix which is used in the lars algorithm.
#' Each signal is weighted 
#'
#' @title lars algorithm for bivariate signal
#' @param CN matrix containing copy-number signals. Each row corresponds to a different signal.
#' @param fracB matrix containing copy-number signals. Each row corresponds to a different signal.
#' @param y vector containing the response associated to each signal
#' @param weightsCN vector of length nrow(CN); weights associated to each signal for the copy-number signal
#' @param weightsFracB vector of length nrow(fracB); weights associated to each signal for the copy-number signal
#' @param meanCN value for centering the copy-number signal (default value = 2)
#' @param maxSteps maximum number of steps for the lars algorithm
#' @param eps tolerance
#' 
#' @return a LarsPath object
#' 
#' @author Quentin Grimonprez
#' 
#' @export
HDlarsbivariate=function(CN,fracB,y,weightsCN=1/apply(CN,1,sd),weightsFracB=1/apply(fracB,1,sd),meanCN=2,maxSteps,eps)
{
  #CN
  if(missing(CN))
    stop("CN is missing.")
  if(!is.numeric(CN) || !is.matrix(CN))
    stop("CN must be a matrix of real.")
  #fracB
  if(missing(fracB))
    stop("fracB is missing")
  if(!is.numeric(fracB) || !is.matrix(fracB))
    stop("fracB must be a matrix of real.")
  #y
  if(missing(y))
    stop("y is missing.")
  #weightsCN
  if(!is.numeric(weightsCN) && !is.vector(weightsCN))
    stop("weightsCN must be a vector of real.")
  if(length(weightsCN)!=nrow(CN))
    stop("The length of weightsCN doesn't match with the number of rows of CN.")
  #weightsFracB
  if(!is.numeric(weightsFracB) && !is.vector(weightsFracB))
    stop("weightsFracB must be a vector of real.")
  if(length(weightsFracB)!=nrow(fracB))
    stop("The length of weightsFracB doesn't match with the number of rows of fracB.")
  #meanCN
  if(!is.double(meanCN))
    stop("meanCN must be a real.")
  
  #y,maxSteps,eps will be checked in HDlars
  
  X=(CN-meanCN)*weightsCN+fracB*weightsFracB
  
  res=HDlars(X,y,maxSteps,eps)
  
  return(res)
}


# setwd("~/Documents//aromaRb")
# folder=c("INCA01-09","INCA10-29","INCA30-39","INCA40-49","INCA50-59","INCA60-69")
# normalTumorArrayINCA=read.csv("normalTumorArrayINCA.csv")
# greffe=which(normalTumorArrayINCA$greffe==1)
# na=which(is.na(normalTumorArrayINCA$greffe))
# etrange=c(2,36)#66 semble greffé #7 CN étrange
# omit=union(union(na,greffe),etrange)
# i=1
# fracB=CN=sampleNames=c()
# for(i in 1:length(folder))
# {
#   CNtemp=getCopyNumberSignal(folder[i],22,normalTumorArrayINCA,onlySNP=TRUE,verbose=FALSE)
#   CN=cbind(CN,CNtemp$copynumber)
#   fracB=cbind(fracB,getFracBSignal(folder[i],22,normalTumorArrayINCA,verbose=FALSE)$fracB$tumor)
#   sampleNames=c(sampleNames,CNtemp$sampleNames)
# }
# position=CNtemp$position
# CN=CN[,-omit]
# fracB=fracB[,-omit]
# fracB=apply(fracB,2,symmetrizeFracB)
# sampleNames=sampleNames[-omit]
# response=normalTumorArrayINCA$rechute[-omit]

HDlarsbivariate2=function(CN,fracB,y,weightsCN=1/apply(CN,1,sd),weightsFracB=1/apply(fracB,1,sd),meanCN=2,maxSteps,eps)
{
  #CN
  if(missing(CN))
    stop("CN is missing.")
  if(!is.numeric(CN) || !is.matrix(CN))
    stop("CN must be a matrix of real.")
  #fracB
  if(missing(fracB))
    stop("fracB is missing")
  if(!is.numeric(fracB) || !is.matrix(fracB))
    stop("fracB must be a matrix of real.")
  #y
  if(missing(y))
    stop("y is missing.")
  #weightsCN
  if(!is.numeric(weightsCN) && !is.vector(weightsCN))
    stop("weightsCN must be a vector of real.")
  if(length(weightsCN)!=nrow(CN))
    stop("The length of weightsCN doesn't match with the number of rows of CN.")
  #weightsFracB
  if(!is.numeric(weightsFracB) && !is.vector(weightsFracB))
    stop("weightsFracB must be a vector of real.")
  if(length(weightsFracB)!=nrow(fracB))
    stop("The length of weightsFracB doesn't match with the number of rows of fracB.")
  #meanCN
  if(!is.double(meanCN))
    stop("meanCN must be a real.")
  
  #y,maxSteps,eps will be checked in HDlars
  
  X=cbind((CN-meanCN)*weightsCN,fracB*weightsFracB)
  
  res=HDlars(X,y,maxSteps,eps)
  
  return(res)
}
