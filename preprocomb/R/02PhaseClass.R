#' @include 01DataClass.R
NULL

## PHASE

setClass("PhaseClass", representation(objectname="character", preprotransformations="list", preimpute="logical"))

#' setphase
#'
#' setphase is a constructor function for defining a phase.
#' setphases initializes a PhaseClass object.
#'
#' @param phasename (character) name of the phase
#' @param preprocessor (character) vector of preprocessors (see ?setpreprocessor) belonging to the phase
#' @param preimpute (logical) whether phase is executed before missing value imputation
#' @return a PhaseClass object
#' @examples
#' ## imputation <- setphase("imputation", c("naomit", "meanimpute"), TRUE)
#' @export
#' @details All elements of argument 'preprocessor' must point to PreprocessorClass objects constructed with function 'setpreprocessor()'.

setphase <- function(phasename, preprocessor, preimpute){

  if (class(phasename)!="character") {stop("Argument 'phasename' must be a character string.")}
  if (class(preprocessor)!="character") {stop("Argument 'preprocessor' must be a character vector.")}
  if (length(preprocessor)==0) {stop("Argument 'preprocessor' must have one or more elements.")}
  if (class(preimpute)!="logical") {stop("Argument 'preimpute' must be a logical (TRUE/FALSE).")}

  listofpreprocessors <- as.list(preprocessor)
  if (any(unlist(lapply(listofpreprocessors, function(x) extends(x, "PreprocessorClass")))=="FALSE")) {
    stop("All elements of argument 'preprocessor' must point to PreprocessorClass objects constructed with function 'setpreprocessor'.") }

  phaseclassobject <- new("PhaseClass", objectname=phasename, preimpute=preimpute)
  phaseclassobject@preprotransformations <- listofpreprocessors
  return(phaseclassobject)
}




