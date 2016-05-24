learnParameters <- function(object) {
  
  # Extract observations and prediction locations
  inputs = object$observations
  pred   = object$predictionLocations
  
  # Put data into an easy parseable format for the backend C++ code
  x = coordinates(inputs)
  
  # JOS: Allow for other names for the dependent variable
  depVar = as.character(object$formulaString[[2]])
  y = as.numeric(inputs@data[[depVar]])
  
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
  
  # Retrieve metadata
  metaData = object$obsChar
  
  # If no metadata, attempt to build it from the variances
  # (observations$oevar) and, if available, biases (observations$oebias)
  # in the observations data structure
  if (is.null(metaData)) 
  {
    metaData <- buildMetadata(inputs);
  }
  
  #------------------------------------------
  # Estimate range/sill using variogram model
  #------------------------------------------
  vario = array()
  vmodels = c("Exp","Gau")   # Restrict to exponential or gaussian models
  varioModel = autofitVariogram(object$formulaString,
                    object$observations,model=vmodels)$var_model
  vario[1] = varioModel$model[2]        # Model used (Gaussian: 1, Exp: 2)
  vario[2] = varioModel[[3]][2]         # Range
  vario[3] = varioModel[[2]][2]         # Sill
  vario[4] = varioModel[[2]][1]         # Nugget
  
  
  #------------------------------------------
  # call the C code
  #------------------------------------------
  try(
    r <- .Call("estimateParams", x, y, vario, obsErrId, sensorId, metaData,
              PACKAGE = "psgp")
  )
  
  #------------------------------------------
  # Retrieve and store PSGP parameters
  #------------------------------------------
  # PSGP covariance parameters (log transformed)
  # We store them in the variogram model - this is a hack really, we just
  # use the variogram model data frame for storage!
  
  # Create dummy variogram model - this allows us to store 2x9 values
  # as the first column is for text
  vmodel = vgm(0, "Exp", 1, 0)
  # Use N/A to indicate that something is "wrong" with this variogram model
  # This should make sure it is not used or returns an error if it is.
  vmodel[1:2,1] = NA    
  vmodel[1,2:9] = r[1:8]           # Store PSGP parameters in columns
  vmodel[2,2:9] = r[9:16]          # 2 to 9
  
  object$variogramModel = vmodel
  object
}

