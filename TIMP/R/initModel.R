"initModel" <-
  # initializes the model class based on a match of the named argument "mod_type" to one of the supported model types
  # Variables directly affected by this function:
  # Initialized: model@clpType
  # Modified: model@clpType, model@datCall
  function (...) 
  {
    dots <- list(...)
    if ("mod_type" %in% names(dots)) {
      if (dots$mod_type == "kin") 
        model <- kin(...)
      if (dots$mod_type == "spec") 
        model <- spec(...)
      if (dots$mod_type == "mass") 
        model <- mass(...)
      if (dots$mod_type == "amp") 
        model <- amp(...)
    }
    if(model@mod_type == "spec")
      model@clpType <- "x"
    else model@clpType <- "x2"
    model@datCall <- append(model@datCall, match.call())   
    model <- initOneModel(model) 
    model
  }