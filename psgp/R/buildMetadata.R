buildMetadata <- function(observations) 
{
  # Retrieve variance for each error model
  obsErrVar  = as.double(observations$oevar)
  
  # If no variances provided, then there is no observation error
  # Otherwise, check for bias
  if (!is.null(obsErrVar)) 
  {
    # Retrieve bias for each error model - this can be omitted by the user
    obsErrBias = as.double(observations$oebias)
    
    # If no bias provided, use zero bias for each error model
    if (is.null(obsErrBias) || length(obsErrBias) == 0) {
      obsErrBias = rep(0.0, length(obsErrVar))
    }
    
    # Check that we have enough variances and biases, i.e.
    # we need at least max(observations$oeid) error models where oeid
    # is a vector of indexes giving the error model used for each
    # observation
    nErrorModels = max(observations$oeid,0)
    if ( length(obsErrVar) < nErrorModels || length(obsErrBias) < nErrorModels )
    {
      print("Observation error index exceeds number of error models")
      print("Make sure that the length of observations$oevar (and observations$oebias ")
      print("if provided) is at least equal to the number of error ")
      print("models requested (as given by max(observations$oeid))")
      print("Observation error is ignored in the rest of the procedure")
    }
    else 
    {
      # Build list of (Gaussian) model descriptors 
      metaData = mapply(FUN = function(...) paste("GAUSSIAN",...,sep=","), 
                        obsErrBias,obsErrVar)
      
      # Append sensor model.
      # The sensor model goes at the end of the metadata and, for the 
      # moment, should always be empty (this is what the C++ library 
      # expects). In the future, when we add support for sensor models,
      # this will allow MathML expressions (if I understand correctly)
      # to be passed in.
      sensorModel = length(obsErrVar) + 1
      metaData[sensorModel] = ""
      
      # Convert vector to list
      metaData = as.list(metaData);
    }
  }
  metaData
}
