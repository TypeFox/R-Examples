
setClass(Class="FDatFitNormal",
         contains="ForecastData",
         representation=representation(
           modelWeights="numeric",
           modelParams="array",
           variance="numeric",
           logLik="numeric",
           exp="numeric",
           tol="numeric",
           maxIter="numeric",
           method="character",
           predType="character",
           iter="numeric",
           model="character",
           modelResults = "list",
           useModelParams = "logical",
           call="call"
           )
         )


