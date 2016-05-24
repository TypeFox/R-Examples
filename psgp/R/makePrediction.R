makePrediction <- function(object, vario)
{
  inputs = object$observations
  pred = object$predictionLocations
  
  # put data into an easy parseable format for the backend C++ code
  x = coordinates(inputs)
  y = as.vector(inputs$value)
  
  ### error variance vector
  e = as.vector(inputs$var)
  
  tx = coordinates(pred)
  
  ### error variance vector
  e = as.vector(inputs$value)
  
  #-------------------------------------
  # Extract observation error components
  #-------------------------------------
  # To each observation corresponds:
  # - the index of an error model (oeid), i.e. the model oeid(i) is used 
  #   for obs i
  # - a sensor model (sensorid) - THIS IS NOT USED YET
  # - and either:
  #   - a full metadata description of the error models, i.e. a list of
  #     strings in the form "<distname>,<bias>,<variance>". For example:
  #     [ "GAUSSIAN,0.0,1.3"
  #       "GAUSSIAN,0.2,1.6" ]
  #     At the moment, only a GAUSSIAN distribution with zero bias is 
  #     allowed by PSGP.
  #   or:
  #   - the variances of the error models (oevar) - which variance is used 
  #     for a particular observation is determined by the index in oeid.
  #   - the biases of the error models (oebias) - same as above for the bias
  #     The variance and bias terms are only taken into account if no metadata
  #     is provided, and are converted to a valid metadata table.
  obsErrId = as.integer(inputs$oeid)
  sensorId = as.integer(inputs$sensor)
  
  # If a metadata has been provided, pass it to PSGP directly
  metaData = object$obsChar
  
  # Otherwise, check if observation error information has been 
  # provided instead
  if (is.null(metaData)) 
  {
    metaData <- buildMetadata(inputs);
  }

  try(
  r <- .Call("predict", x, y, tx, vario, obsErrId,
             sensorId, metaData, PACKAGE = "psgp")
  )
}


