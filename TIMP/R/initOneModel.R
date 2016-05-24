"initOneModel" <-
  # Initializes some fields of the model class, calls the generic method initModelClass which depends on mod_type
  function (model) 
  {
    model@weight <- length(model@weightpar) != 0
    model@lclp0 <- length(model@clp0) != 0
    model@lclpequ <- length(model@clpequspec) != 0
    model <- initModelClass(model)    
    if(length(model@clp0) != 0) 
      model@usecompnames0 <- is.character(unlist(lapply(model@clp0, 
                                                        function(x){x[["comp"]]} ))) 
    if(length(model@clpequspec) != 0) 
      model@usecompnamesequ <- is.character(unlist(lapply(model@clpequspec, 
                                                          function(x){x[["to"]]} ))) 
    if(model@usecompnames0 || model@usecompnamesequ) {
      if(length(model@compnames) < model@ncomp) {
        model@compnames <- append(model@compnames, 
                                  1:(model@ncomp-length(model@compnames)))
      }
      if(length(model@compnames > model@ncomp))
        cat(paste("WARNING: too many component names, last", 
                  model@compnames-model@ncomp, "names not used\n"))
    }
    if(length(model@autoclp0) > 0){
      model@clp0 <- getAutoClp0(model)
      model@autoclp0 <- list() ## don't want to store the old result
    }
    if(model@mod_type=="kin") {
      if(model@irffun=="doublegaus") {
        model@doublegaus <- TRUE
      }
      if(model@irffun=="multiplegaus") {
        model@multiplegaus <- TRUE
      }
      if(model@dispmu && length(model@parmu[[1]])>0 && length(model@lambdac)==0){
        warning(sprintf("model@dispmu = %s but lamdac not specified!",model@dispmu))
      }
    }
    model@clpCon <- getClpConstr(
      if (model@clpType == "x")
        clp = model@x else clp = model@x2,
      clp0 = model@clp0, 
      clpequspec = model@clpequspec, 
      ncomp = model@ncomp, 
      compnames = model@compnames)
    model@x <- model@x * model@scalx
    model@fvecind <- getFixed(model)
    model@pvecind <- getPrel(model)
    model@mvecind <- model@nvecind <- getMvec(model)	
    model <- getConstrained(model)
    model <- addPrel(model)
    model
  }
