## ###################################################
## Function to stratifiy data from a data.frame-object
## ###################################################
ps.makestrata.data.frame <- function(object,                 
                                     breaks             = NULL,            
                                     name.stratum.index = "stratum.index",   
                                     stratified.by      = NULL,  
                                     ...)
{
  ## ############
  ## Check object
  if (missing(object)){
    stop("Argument 'object' is missing.")
  }else{
    if(!inherits(object,"data.frame")){
      stop("Argument 'object' is not of class 'data.frame'.")
    }else{
      data <- object
    }
  }
  
  ## ###################
  ## Check stratified.by
  if ( is.null(stratified.by) ){
    stop("Argument 'stratified.by' is needed.")
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
  if(any(names(data) == name.stratum.index))
    stop(paste("Argument 'name.stratum.index'=",
               name.stratum.index, " already exists in data.", sep=""))
  
  ## ############################
  ## Check breaks and name.pscore
  if(!is.null(breaks))
    if (length(breaks)!=1)
      if (any(strata.vec > max(breaks)) | any(strata.vec < min(breaks)))
        warning("Either any(data[,stratified.by] > max(breaks)) or any(data[,stratified.by] < min(breaks)) holds. NA values in strata results!")
  
  ## ###########
  ## Make strata
  if(!is.null(breaks)){
    if(length(breaks)!=length(unique(breaks))){
      stop(paste("Argument 'breaks' =", breaks,
                 "are not unique", sep=""))
    }else{
      strata=cut(strata.vec, breaks, incl=TRUE,...)
    }
  }else{
    strata <- factor(round(strata.vec,3))
  }

  intervals <- levels(strata)
  levels(strata) <- c(1:length(intervals))
 
  data[,name.stratum.index] <- strata
  
  stra.df <- list(data               = data,
                  intervals          = intervals,
                  stratum.index      = strata,
                  name.stratum.index = name.stratum.index,
                  stratified.by      = strata.name)
  
  class(stra.df) <- c("stratified.data.frame")

  return(stra.df)
}
