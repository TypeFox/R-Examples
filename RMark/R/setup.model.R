#' Defines model specific parameters (internal use)
#' 
#' Compares \code{model}, the name of the type of model (eg "CJS") to the list
#' of acceptable models to determine if it is supported and then creates some
#' global fields specific to that type of model that are used to modify the
#' operation of the code.
#' 
#' In general, the structure of the different types of models (e.g.,
#' "CJS","Recovery",...etc) are very similar with some minor exceptions.  This
#' function is not intended to be called directly by the user but it is
#' documented to enable other models to be added.  This function is called by
#' other functions to validate and setup model specific parameters.  For
#' example, for live/dead models, the length of the capture history is twice
#' the number of capture occasions and the number of time intervals equals the
#' number of capture occasions because the final interval is included with dead
#' recoveries.  Whereas, for recapture models, the length of the capture
#' history is the number of capture occasions and the number of time intervals
#' is 1 less than the number of occasions.  This function validates that the
#' model is valid and sets up some parameters specific to the model that are
#' used in the code.
#' 
#' @param model name of model type (must be in vector \code{valid.models})
#' @param nocc length of capture history string
#' @param mixtures number of mixtures
#' @export
#' @return model.list - a list with following elements \item{etype}{encounter
#' type string for MARK input; typically same as model} \item{nocc}{number of
#' capture occasions} \item{num}{number of time intervals relative to number of
#' occasions (0 or -1)} \item{mixtures}{number of mixtures if any}
#' \item{derived}{logical; TRUE if model produces derived estimates}
#' @author Jeff Laake
#' @seealso \code{\link{setup.parameters}}, \code{\link{valid.parameters}}
#' @keywords utility
setup.model <-
function(model,nocc,mixtures=1)
{
#
# Value: 
#
#   model.list - a list with following elements
#                  etype - encounter type string; typically same as model name
#                  nocc  - number of capture occasions
#                  num   - number of time intervals relative to number of occasions (0 or -1)
#                  mixtures - number of mixtures if any
#                  derived - TRUE if model produces derived parameters
#
#
# Read in parameter definitions
	fdir=system.file(package="RMark")	
	fdir=file.path(fdir,"models.txt")	
	model_definitions=read.delim(fdir,header=TRUE,
			colClasses=c("numeric","character","character",rep("logical",4),rep("numeric",3),"logical"))
    model_def=model_definitions[model_definitions$model==model,]	
    if(nrow(model_def)==0)
        stop("Invalid type of model = ",model," Valid types are\n", paste(model_definitions$model,collapse="\n"))
	if(mixtures==1) 
		model_def$mixtures=model_def$default.mixtures
	else
		model_def$mixtures=mixtures
	model_def$default.mixtures=NULL
	model_def$nocc=nocc/model_def$divisor
	if(model_def$derived)
	{
		fdir=system.file(package="RMark")	
		fdir=file.path(fdir,"DerivedPar.txt")	
		deriv_pars=read.delim(fdir,header=TRUE,	colClasses=c("numeric","character"))
		model_def$derived_labels=list(deriv_pars$dpar_label[deriv_pars$MarkNumber==model_def$MarkNumber])
	}
    return(as.list(model_def))
}
