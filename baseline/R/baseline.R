## $Id: baseline.R 192 2012-06-19 08:36:53Z kristl $
### Main baseline correction function, and definition of class baseline.

###
### Baseline class
###
setClass("baseline",
         representation(baseline = "matrix", corrected = "matrix",
                        spectra = "matrix", call = "language")
         )

###
### Top level baseline correction function
###

baseline <- function (spectra, method = "irls", ...) {
    ## Get baseline algorithm function name:
	if(exists("baselineAlgorithms",envir=.GlobalEnv)){
		bA <- get("baselineAlgorithms",envir=.GlobalEnv)
	} else {
		bA <- baselineAlgorithms
	}
    method <- match.arg(method, names(bA))
    baseFunc <- funcName(bA[[method]])

    ## Run baseline algorithm:
    res <- do.call(baseFunc, list(spectra, ...))

    ## Build and return the object:
    new("baseline",
        baseline = res$baseline,
        corrected = res$corrected,
        spectra = spectra,
        call = match.call()
        )
}


###
### Extraction methods
###
setGeneric("getSpectra", function(object) standardGeneric("getSpectra"))
setMethod("getSpectra", "baseline", function(object) object@spectra)
setGeneric("getCorrected", function(object) standardGeneric("getCorrected"))
setMethod("getCorrected", "baseline", function(object) object@corrected)
setGeneric("getBaseline", function(object) standardGeneric("getBaseline"))
setMethod("getBaseline", "baseline", function(object) object@baseline)
setGeneric("getCall", function(object) standardGeneric("getCall"))
setMethod("getCall", "baseline", function(object) object@call)
