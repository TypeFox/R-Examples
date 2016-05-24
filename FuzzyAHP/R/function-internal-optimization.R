######################################################
# this file contains internal functions
# for calculation of fuzzy weights
######################################################

setGeneric(".weightsLimits",
           function(data, i) standardGeneric(".weightsLimits"))

setMethod(
  f=".weightsLimits",
  signature(data = "FuzzyPairwiseComparisonMatrix", i = "numeric"),
  definition=function(data, i)
  {
    p = nrow(data@fnMin)

    # calculate lower value of weight
    prodL = prod(data@fnMin[i,])^(1/p)

    optL1 = .optProLower(data, i, "min")
    optL2 = .optProLower(data, i, "max")

    wL = prodL /(prodL + max(optL1, optL2))

    # calculate upper value of weight
    prodU = prod(data@fnMax[i,])^(1/p)

    optU1 = .optProUpper(data, i, "min")
    optU2 = .optProUpper(data, i, "max")

    wU = prodU /(prodU + min(optU1, optU2))

    return(c(wL,wU))
  }
)

# function for calculating upper limits
setGeneric(".optProUpper",
           function(data, i, type) standardGeneric(".optProUpper"))

setMethod(
  f=".optProUpper",
  signature(data = "FuzzyPairwiseComparisonMatrix", i = "numeric", type = "character"),
  definition=function(data, i, type)
  {
    sum = 0
    p = nrow(data@fnMin)

    element = data@fnMin

    if(type == "min"){
      matrix = data@fnMin

    }
    else if(type == "max"){
      matrix = data@fnMax
    }
    else{
      stop(paste("Unrecognized type (should be min or max) is: ",type,".", sep = ""))
    }

    return(.optPro(data, i, matrix, element))
  }
)

# function for calculating lower limits
setGeneric(".optProLower",
           function(data, i, type) standardGeneric(".optProLower"))

setMethod(
  f=".optProLower",
  signature(data = "FuzzyPairwiseComparisonMatrix", i = "numeric", type = "character"),
  definition=function(data, i, type)
  {
    sum = 0
    p = nrow(data@fnMin)

    element = data@fnMax

    if(type == "min"){
      matrix = data@fnMin

    }
    else if(type == "max"){
      matrix = data@fnMax
    }
    else{
      stop(paste("Unrecognized type (should be min or max) is: ",type,".", sep = ""))
    }

    return(.optPro(data, i, matrix, element))
  }
)

# optimization function
setGeneric(".optPro",
           function(data, i, matrix, element) standardGeneric(".optPro"))

setMethod(
  f=".optPro",
  signature(data = "FuzzyPairwiseComparisonMatrix", i = "numeric", matrix = "matrix", element = "matrix"),
  definition=function(data, i, matrix, element)
  {
    sum = 0
    p = nrow(matrix)

    for (k in 1:p){
      if(k == i){
        next
      }

      if(k>1){
        prod1 = 1
        for(l in 1:(k-1)){
          if(l == i){
            next
          }
          prod1 = prod1 * (1/matrix[l,k]) #prod(1/matrix[l,k]) # remove i from set
        }
      }else{
        prod1 = 1
      }

      if((k+1)<=p){
        prod2 = 1
        for(l in (k+1):p){
          if(l == i){
            next
          }
          prod2 = prod2 * matrix[k,l] #prod(matrix[k,(k+1):p]) # remove i from set
        }
      }else{
        prod2 = 1
      }

      sum = sum + (element[k,i] * prod1 * prod2)^(1/p)
    }

    return(sum)
  }
)
