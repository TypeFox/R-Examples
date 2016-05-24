#' LHS Accessor Functions.
#' 
#' 
#'	Instead of using the $ operator, using these accessor functions 
#'	is the preferred method for accessing the data and result
#'	data frames from an \code{\link{LHS}} or \code{\link{PLUE}} object, 
#'  as the internal structure of
#'	the object may vary between versions of the package.
#'
#'	\code{get.data} returns a data.frame consisting on the input data.
#'
#'	\code{get.results} returns an array with the model results. See the 
#'	vignette on multiple runs for details on the \code{get.mean} argument.
#'
#'	\code{get.N}, \code{get.ninputs}, \code{get.noutputs} return a single number each,
#'	with the number of points in the hypercube (or sampling), number of input factors and number of
#'	response variables.
#'
#'	\code{get.repetitions} returns the number of model repetitions for each data point,
#'	created by \code{LHS(model, factors, N, repetitions=X)}, or by \code{tell}ing several
#'	result sets to the same LHS object (or PLUE object).
#'
#' @param obj The LHS or PLUE object
#' @param get.mean In case of stochastic models, when several model runs are required for the
#'	  same data point, the \code{data} slot of the LHS object contains all the 
#'	  model outputs. Use \code{get.mean=TRUE} to get the average values for each point,
#'	  or \code{get.mean=FALSE} to get all the results.

#' @rdname accessors
#' @export
get.results <-
function(obj, get.mean=TRUE) { 
	# IF the object is incomplete (see tell method)
	if (is.null(dim(obj$res))) {return (NA)}

	if(get.mean) {
		return (apply(obj$res, c(1,2), mean)) 
	} else {
		return (obj$res) 
	}
}

#' @rdname accessors
#' @export
get.data <-
function(obj) {
	return (obj$data) 
}

#' @rdname accessors
#' @export
get.N <- function(obj) {return (dim(get.data(obj))[1]) }

#' @rdname accessors
#' @export
get.ninputs <- function(obj) {return (dim(get.data(obj))[2]) }

#' @rdname accessors
#' @export
get.noutputs <- function(obj) {return (dim(get.results(obj))[2]) }

#' @rdname accessors
#' @export
get.repetitions <- function(obj) {return (dim(get.results(obj, get.mean=FALSE))[3]) }
