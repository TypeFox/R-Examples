if (!isGeneric("createDataPartition")) {
  setGeneric("createDataPartition")
}

if (!isGeneric("createResample")) {
  setGeneric("createResample")
}

if (!isGeneric("createFolds")) {
  setGeneric("createFolds")
}

if (!isGeneric("createMultiFolds")) {
  setGeneric("createMultiFolds")
}

setMethod("createDataPartition", signature(y = ".CaretHyperspectral"),
          definition = function(y,
                                times = 1,
                                p = 0.5,
                                list = TRUE,
                                groups = min(5, length(y)))
{  
  return(createDataPartition(.getResponseVar(y, 
                                              advice = c("createDataPartition", "setResponse")),
                             times = times,
                             p = p,
                             list = list,
                             groups = groups))
  
})

setMethod("createResample", signature(y = ".CaretHyperspectral"),
          definition = function(y,
                                times = 10, 
                                list = TRUE)
{  
  return(createResample(.getResponseVar(y, 
                                         advice = c("createResample", "setResponse")),
                        times = times,
                        list = list
                       ))
  
})

setMethod("createFolds", signature(y = ".CaretHyperspectral"),
          definition = function(y,
                                k = 10, 
                                list = TRUE, 
                                returnTrain = FALSE)
{  
  return(createFolds(.getResponseVar(y, 
                                      advice = c("createFolds", "setResponse")),
                     k = k,
                     list = list,
                     returnTrain = returnTrain
                    ))
  
})

setMethod("createMultiFolds", signature(y = ".CaretHyperspectral"),
          definition = function(y,
                                k = 10,
                                times = 5)
{  
  return(createMultiFolds(.getResponseVar(y, 
                                           advice = c("createMultiFolds", "setResponse")),
                          times = times,
                          k = k
                         ))
  
})


