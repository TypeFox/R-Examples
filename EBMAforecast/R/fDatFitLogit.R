
setClass(Class="FDatFitLogit",
         contains="ForecastData",
         representation=representation(
           modelWeights="numeric",
           modelParams="array",
           logLik="numeric",
           exp="numeric",
           tol="numeric",
           maxIter="numeric",
           method="character",
           iter="numeric",
           model="character",
           modelResults = "list",
           useModelParams = "logical",
           call="call"
           ),
         validity=function(object){
           if(length(object@modelWeights)>0){
             if(sum(object@modelWeights)<=.99 || sum(object@modelWeights)>1.01){
             stop("Model weights should sum to approximately one.")
           }
           }
         }
         )

