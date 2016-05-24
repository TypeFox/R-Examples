MgCheck <- function(object,model.args,model.parms,SplitMethod,verbose=TRUE){
    checkF <- lapply(1:length(object),function(f){
        fit <- object[[f]]
        if(SplitMethod != "noinf" && is.null(fit$call))
            stop(paste("Element",f," '",names(object)[f],"' does not have a call argument\nThis is needed to rebuild the model during cross-validation."))
        else fit$call$data <- NULL
    })
  
  # check model.args
  # --------------------------------------------------------------------
  if (!is.null(model.args)){
    if (!(is.list(model.args))){
      warning(paste("Argument model.args is not a list and therefore ignored." ))
      model.args <- NULL
    }
    else{
      if (!(all(match(make.names(names(model.args),unique=TRUE),names(object),nomatch=FALSE)))){
        if (verbose==TRUE)
          warning(paste("model.args should be a named list matching the entries of the object. Assume now that they are given in the correct order" ))
      }
      else{
        model.args <- model.args[names(object)]
      }
    }
  }
    
  # check model.parms
  # --------------------------------------------------------------------
  if (!is.null(model.parms)){
    if (!(is.list(model.parms))){
      warning(paste("Argument model.parms is not a list and therefore ignored." ))
      model.args <- NULL
    }
    else{
      if (!(all(match(make.names(names(model.parms),unique=TRUE),names(object),nomatch=FALSE)))){
        if (verbose==TRUE)
          warning(paste("model.parms should be a named list matching the list of model.\nIt is assumed that they are given in the correct order" ))
      }
      else{
        model.parms <- model.parms[names(object)]
      }
    }
  }
  list(model.args=model.args,model.parms=model.parms)
}
