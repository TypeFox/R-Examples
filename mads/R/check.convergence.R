#' Checks whether the model has converged 
#'  
#' @param ddf.model ddf object
#' @return boolean
#' @author Laura Marshall
#'
check.convergence <- function(ddf.model){
  model.type <- ddf.model$method
  converged <- NULL
  if(model.type%in%c("ds")){
    if(ddf.model$ds$converge == 0){
      converged <- TRUE
    }else{
      converged <- FALSE
    }
  }else if(model.type%in%c("trial", "trial.fi", "io", "io.fi")){
    converged <- ddf.model$mr$mr$converged
  }
  return(converged)
}