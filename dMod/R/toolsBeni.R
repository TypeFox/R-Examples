
#' Search for symmetries in the loaded model
#' 
#' @param f eqnList object or vector containing ODEs
#' @param obsvect vector of observation functions
#' @param prediction vector containing prediction to be tested
#' @param initial vector containing initial values
#' @param ansatz type of infinitesimal ansatz used for the analysis (uni, par, multi)
#' @param pMax maximal degree of infinitesimal ansatz
#' @param inputs specify the input variables
#' @param fixed variables to concider fixed
#' @param cores maximal number of cores used for the analysis
#' @param allTrafos do not remove transformations with a common parameter factor
#' @return NULL
#' @export
symmetryDetection <- function(f, obsvect = NULL, prediction = NULL, initial = NULL, ansatz = 'uni', pMax = 2, inputs = c(), fixed = c(), cores = 1, allTrafos = FALSE){
  
  if (is.element("observables", names(attributes(f)))){
    f <- f[1:(length(f)-length(attr(f,"observables")))]
  }
  f <- as.character(lapply(1:length(f), function(i) paste(names(f)[i],'=',f[i])))
  
  obsvect <- as.character(lapply(1:length(obsvect), function(i) paste(names(obsvect)[i],'=',obsvect[i])))
  prediction <- as.character(lapply(1:length(prediction), function(i) paste(names(prediction)[i],'=',prediction[i])))
  initial <- as.character(lapply(1:length(initial), function(i) paste(names(initial)[i],'=',initial[i]))) 
  
  rPython::python.load(paste(system.file(package="R2CdeSolve"),"/code/functions.py", sep = ""))
  rPython::python.load(paste(system.file(package="R2CdeSolve"),"/code/readData.py", sep = ""))
  rPython::python.load(paste(system.file(package="R2CdeSolve"),"/code/buildSystem.py", sep = ""))
  rPython::python.load(paste(system.file(package="R2CdeSolve"),"/code/symmetryDetection.py", sep = ""))
  
  rPython::python.call("symmetryDetection", f, obsvect, prediction, initial, ansatz, pMax, inputs, fixed, cores, allTrafos)
}

#' Do a variable transformation in the ODE
#' 
#' @param observables Named character vector. The names are the new variable names, the vector
#' entries define the new variables in terms of the old ones.
#' @param f An object of class \code{eqnList}, see \link{generateEquations}.
#' @param dynvar Character vector with the old variable names
#' @param stoi The stoichiometric matrix
#' @param flows Character vector with the rate expressions
#' @param conserved Logical. If true, the conserved quantities derived from the
#' stoichiometric matrix are automatically used for extending the vector of observables. See details.
#' @details Usually, the function is called by either using the \code{f} argument and leaving the
#' other arguments \code{NULL} or by leaving \code{f} NULL and defining the ODE by the arguments
#' \code{dynvar}, \code{stoi} and \code{flows}.
#' The \code{observables} vector can have less entries than the vector \code{dynvar}. In this
#' case, the observables are automatically extended by old variables to generate a full rank
#' variable transformation. If \code{conserved} is \code{TRUE}, the conserved quantities are
#' preferentially used to extend the observables vector. Consequently, the transformed equations
#' will return a certain number of zero-equation.
#' @return Named character vector with the ODE expressed in the new variables. In addition,
#' attributes "variables" (the variable transformation) and "inverse" (the inverse transformation)
#' are returned.
#' @export
variableTransformation <- function(observables, f = NULL, dynvar = NULL, stoi = NULL, flows = NULL, conserved=TRUE){
  
  observation <- paste(names(observables), "=", observables) 
  if(is.null(dynvar)) dynvar <- attr(f, "species")
  if(is.null(stoi)) stoi <- c(t(attr(f, "SMatrix"))); stoi[is.na(stoi)] <- 0
  if(is.null(flows)) flows <- attr(f, "rates")
  
  
  rPython::python.load(paste(system.file(package="R2CdeSolve"),"/code/functions_obs.py", sep = ""))
  rPython::python.load(paste(system.file(package="R2CdeSolve"),"/code/extendObservation.py", sep = ""))
  
  
  out <- rPython::python.call("getObservation", observation, dynvar, stoi, flows, conserved)
  
  variables <- out[[3]]; names(variables) <- out[[1]]
  f <- out[[4]]; names(f) <- out[[1]]
  inverse <- out[[5]]; names(inverse) <- dynvar
  
  attr(f, "variables") <- variables
  attr(f, "inverse") <- inverse
  
  return(f)
  
}

