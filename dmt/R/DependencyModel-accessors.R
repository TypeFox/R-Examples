# (C) 2008-2011 Leo Lahti and Olli-Pekka Huovilainen
# All rights reserved.
# FreeBSD License (keep this notice)


setMethod("getW","DependencyModel", 
  function(model) { 
    return(model@W) 
  } 
) 

setMethod("getPhi","DependencyModel", 
  function(model) { 
    return(model@phi) 
  } 
) 

setMethod("getScore","DependencyModel", 
  function(model) { 
    return(model@score) 
  } 
) 

setMethod("getParams","DependencyModel", 
  function(model) { 
    return(model@params) 
  } 
) 

setMethod("getModelMethod","DependencyModel", 
  function(model) { 
    return(model@method) 
  } 
) 

setMethod("getWindowSize","DependencyModel", 
  function(model) { 
    if (is.null(getW(model)$X)) {
      return(dim(getW(model))[1])
    }
    else {
      return(c(dim(getW(model)$X)[1],dim(getW(model)$Y)[1]))
    }
  } 
) 

setMethod("getZ","DependencyModel",
  function(model,X,Y) {
    if (length(model@z) > 0) {
      return(model@z)
    } else {
      # Calculate latent variable if it hasn't been calculate with the model    
      if (length(model@data) == 0) {
        # Original data not included with the model
        if (is.null(model@W$Y)) {
          # z for dep model with one data
          # Give error if original data is not given
          if (missing(X)) stop("Original data set is needed as a parameter because latent variable was not calculated while calculating dependency model and original data was not included with the model")
          return(z.expectation(model,X))
        } else {
          # z for dep model with two datas
          # Give error if original data is not given
          if (missing(X) || missing(Y)) stop("Original data sets are needed as parameters because latent variable was not calculated while calculating dependency model and original data was not included with the model")
          return(z.expectation(model,X,Y))     
        }
      } else {
        # Calculate z with included data        
        # z for dep model with one data
        if (is.null(model@W$Y)) return(z.expectation(model,model@data$X))
        # z for dep model with two datas
        else return(z.expectation(model,model@data$X,model@data$Y))
      }
    }
  }
)