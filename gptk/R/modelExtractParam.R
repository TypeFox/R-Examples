modelExtractParam <-
function (model, only.values=TRUE, untransformed.values=FALSE) {
#   if (any(.packages(all.available=TRUE)=="tigre") && is.GPModel(model))
#     model <- modelStruct(model)
  
  funcName <- paste(model$type, "ExtractParam", sep="")
  func <- get(funcName, mode="function")
  params <- func(model, only.values=only.values, untransformed.values=untransformed.values)

  if ( !only.values ) {
    origNames <- names(params)
    if ( "paramGroups" %in% names(model) ) {
      paramGroups <- model$paramGroups
      for ( i in seq(length.out=dim(paramGroups)[2]) ) {
        ind <- grep(1, paramGroups[,i])
        if ( is.list(params) ) {
          names(params)[i] <- origNames[ind[1]]
          for ( j in seq(2, length.out=length(ind)-1) )
            names(params)[i] <- paste(names(params)[i], origNames[ind[j]],sep="/")
        }

        paramGroups[ind[seq(2,length(ind),length=length(ind)-1)], i] <- 0
      }

      params <- params%*%paramGroups
    }
  }
  return (params)
}
