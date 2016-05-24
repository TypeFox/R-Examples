## ###############################################
## Function to stratifiy data from a pscore-object
## ###############################################
ps.makestrata.pscore <- function(object,                
                                 breaks             = NULL, 
                                 name.stratum.index = "stratum.index",  
                                 stratified.by      = NULL,      
                                 ...)
{
  ## #############
  ## Check objects
  if (missing(object)){
    stop("Argument 'object' is missing.")
  }else{
    if(!inherits(object,"pscore")){
      stop("Argument 'object' is not of class 'pscore'.")
    }else{
      data <- object$data
    }
  }

  ## ###################
  ## Check stratified.by
  if ( is.null(stratified.by) ){
    strata.vec  <- object$pscore
    strata.name <- object$name.pscore
  }else{
    if (is.character(stratified.by) | is.numeric(stratified.by)){
      A <- find.sel(data     = data,
                    sel      = stratified.by,
                    sel.name = "stratified.by")
      strata.vec <- A[,1]
      strata.name <- names(A)[1]
    }else{
      stop("Argument 'stratified.by' must be numeric or a string.")
    }
  }

  ## #################
  ## Check name.strata
  if(any(names(data) == name.stratum.index)){
    stop(paste("Argument 'name.stratum.index'=",
               name.stratum.index,
               " already exists in data.", sep=""))
  }

  ## ############################
  ## Check breaks and name.strata
  if (!is.null(breaks))
    if (length(breaks)!=1)
      if (any(strata.vec > max(breaks)) | any(strata.vec < min(breaks)))
        warning("Either any(data[,stratified.by] > max(breaks)) or any(data[,stratified.by] < min(breaks)) holds. NA values in strata results!")
  
  ## ###########
  ## Make strata
  if(!is.null(breaks)){
    if(length(breaks) != length(unique(breaks))){
      stop(paste("Argument 'breaks' =", breaks, "are not unique.", sep=""))
    }else{
      strata <- cut(strata.vec, breaks, incl=TRUE,...)
    }
  }else{
    strata <- factor(round(strata.vec,3))
  }

  intervals <- levels(strata)
  levels(strata) <- c(1:length(intervals))


  ## Output
  object$data[, name.stratum.index] <- strata
  object$stratum.index              <- strata
  object$name.stratum.index         <- name.stratum.index
  object$intervals                  <- intervals
  object$stratified.by              <- strata.name
  
  ##class(object) <- c("stratified.pscore")

  class(object) <- c("stratified.pscore",
                     class(object)[class(object)!="stratified.pscore"])

 
  return(object)
  
}
