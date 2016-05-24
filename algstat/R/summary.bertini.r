#' Summarize Bertini Output
#'
#' This function summarizes the output from Bertini.
#' 
#' @param object an object of class bertini
#' @param ... additional parameters
#' @usage \method{summary}{bertini}(object, ...)
#' @return Invisible string of the printed object.
#' @export
#' @examples
#' 
#' # see ?bertini
#' 
#' 
#' 
summary.bertini <- function(object, ...){
return()	
  ## argument checking and basic variable setting
  stopifnot(is.bertini(object))  
  
  
  ## initialize out list
  out <- as.list(vector(length = 4))
  names(out) <- c("real_finite_solutions", "finite_solutions", 
    "nonsingular_solutions", "singular_solutions")
  
   
  ## find out variables lis
  vars <- str_replace(object$main_data[2], "Variables:  ", "")
  vars <- str_split(vars, " ")[[1]]
  p <- length(vars)
  
  
  ## finite_solutions
  if(length(object$finite_solutions) == 1 && object$finite_solutions == ""){
    out$finite_solutions <- NA
  } else {
  	fSolns <- object$finite_solutions
    fSolns <- fSolns[-c(1,2)]
    fSolns <- fSolns[-length(fSolns)]
    fSolns <- zapsmall(as.numeric(
      sapply(strsplit(fSolns, " "), function(x) x[1])
    ))
    fSolns <- matrix(fSolns[!is.na(fSolns)], ncol = p, byrow = TRUE)
    fSolns <- as.data.frame(fSolns)
    names(fSolns) <- vars
    out$finite_solutions <- fSolns
  } 
  
  
  ## real_finite_solutions
  if(length(object$real_finite_solutions) == 1 && object$real_finite_solutions == ""){
    out$real_finite_solutions <- NA
  } else {
  	rSolns <- object$real_finite_solutions
    rSolns <- rSolns[-c(1,2)]
    rSolns <- rSolns[-length(rSolns)]
    rSolns <- zapsmall(as.numeric(
      sapply(strsplit(rSolns, " "), function(x) x[1])
    ))
    rSolns <- matrix(rSolns[!is.na(rSolns)], ncol = p, byrow = TRUE)
    rSolns <- as.data.frame(rSolns)
    names(rSolns) <- vars
    out$real_finite_solutions <- rSolns
  }  
  
  ## nonsingular_solutions
  if(str_sub(object$nonsingular_solutions[1], 1, 1) == "0"){
    out$nonsingular_solutions <- NA
  } else {
  	nsSolns <- object$nonsingular_solutions
    nsSolns <- nsSolns[-c(1,2)]
    nsSolns <- nsSolns[-length(nsSolns)]
    nsSolns <- zapsmall(as.numeric(
      sapply(strsplit(nsSolns, " "), function(x) x[1])
    ))
    nsSolns <- matrix(nsSolns[!is.na(nsSolns)], ncol = p, byrow = TRUE)
    nsSolns <- as.data.frame(nsSolns)
    names(nsSolns) <- vars
    out$nonsingular_solutions <- nsSolns
  } 
  
  ## singular_solutions
  if(str_sub(object$singular_solutions[1], 1, 1) == "0"){
    out$singular_solutions <- NA
  } else {
  	sSolns <- object$singular_solutions
    sSolns <- sSolns[-c(1,2)]
    sSolns <- sSolns[-length(sSolns)]
    sSolns <- zapsmall(as.numeric(
      sapply(strsplit(sSolns, " "), function(x) x[1])
    ))
    sSolns <- matrix(sSolns[!is.na(sSolns)], ncol = p, byrow = TRUE)
    sSolns <- as.data.frame(sSolns)
    names(sSolns) <- vars
    out$singular_solutions <- sSolns
  }      
  
  out 
    
}



