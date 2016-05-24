if (!isGeneric("preProcess")) {
  setGeneric("preProcess")
}

if (!isGeneric("predict")) {
  setGeneric("predict")
}

setMethod("preProcess", signature(x = ".CaretHyperspectral"),
          definition = function(x, ...)
{
  x <- .getAllPredictors(x, NA)
  x <- x[,sapply(1:ncol(x), 
                 function(i, data) is.numeric(data[,i]), x)]
  x <- preProcess(x, ...)
  return(new(".preProcessHyperspectral", preProcess = x))
})

setMethod("predict", signature(object = ".preProcessHyperspectral"),
          definition = function(object, newdata, ...)
{
  object <- object@preProcess  
  if (!(class(newdata) %in% .getCaretCompatibleClasses()))
    return(predict(object, newdata = newdata, ...))
    
  backup <- newdata  
  newdata_all <- .getAllPredictors(newdata, NA)
  is_num <- sapply(1:ncol(newdata_all),
                   function(i, data) is.numeric(data[,i]), newdata_all)
  newdata <- newdata_all[,is_num]
  
  predicted <- predict(object, newdata = newdata, ...)
  
  newdata_all[,is_num] <- predicted
  
  predicted_spectral <- newdata_all[,attr(newdata_all, "spectral")]
  if (attr(newdata_all, "useattributes"))
    attribute(backup)[,.getCaretParameter(object, "predictor")] <- newdata_all[,attr(newdata_all, "spectral")*(-1)]
  
  if (class(backup) == "Speclib")
    spectra(backup) <- as.matrix(predicted_spectral)
  
  if (class(backup) == "Nri") ### Konvertierung nicht optimal!
  {    
    values <- numeric(length = length(backup@nri@values))
    values[] <- NA
    incr <- length(backup@nri@values)/nrow(predicted_spectral)
    for (i in 1:ncol(predicted_spectral))
    { 
      index <- seq(i, length(values), incr)
      values[index] <- predicted_spectral[,i]
    }
    backup@nri <- distMat3D(values, ncol = ncol(backup@nri), nlyr = nrow(predicted_spectral))
  }
  return(backup)
})

