#' Simplified Optimizer for paleotree Likelihood Functions
#'
#' This function is a deliberately simplistic automater for the function
#' \code{\link{optim}} and the use of the "L-BFGS-B" optimizing method,
#' with initial parameter values and bounds provided with \code{\link{parInit}},
#' \code{\link{parLower}} and \code{\link{parUpper}}. It is mainly provided here
#' as a shorthand to be used in educational demonstrations where model-fitting
#' is not the primary focus, and use in actual analyses should be avoided.
#'
#' @details
#' This is mainly provided in this publicly released package for pedagogical
#' reasons. Users seeking an optimizer for their own analytical purposes
#' should write their own optim function.

#' @param modelFun A likelihood function for a model, of class 'paleotreeFunc'.

#' @return
#' Returns the results from using optim.

#' @seealso
#' \code{\link{constrainParPaleo}} and \code{\link{modelMethods}}

#' @examples
#' 
#' # This function simply replicates optim() as shown below
#'     # where modelFun is the likelihood function
#' 
#' #optim(parInit(modelFun),modelFun,
#' #		lower=parLower(modelFun),upper=parUpper(modelFun), 
#' #		method="L-BFGS-B",control=list(maxit=1000000))

#' @export
optimPaleo<-function(modelFun){
	if(!is(modelFun,'paleotreeFunc')){
		stop("Given function does not appear to be a paleotree likelihood function")}
	#
	res<-optim(parInit(modelFun),modelFun,
		lower=parLower(modelFun),upper=parUpper(modelFun), 
		method="L-BFGS-B",control=list(maxit=1000000))
	return(res)
	}