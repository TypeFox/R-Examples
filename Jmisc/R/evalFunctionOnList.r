#' This function evaluates a function \code{x} under an environment which is created by a list. All elements of the list is local to the function; other words all elements of the list can be accessed directly by the function. 
#' A new environment is created and each element of \code{variables} is assigned to the new environment. Then the environment associated with the \code{x} is updated with the new environment. Finally \code{x(...)} is evaluated and return the result.
#' @name evalFunctionOnList
#' @aliases evalFunctionOnList
#' @title Evaluate Function Under Local Variables
#' @param x A function to be called
#' @param variables A list to be converted to an environment
#' @param \dots Further arguments to \code{x}
#' @param parent_env parent environment
#' @return Return value of the \code{x(\dots)}.
#' @author TszKin Julian Chan \email{ctszkin@@gmail.com}
#' @seealso \code{\link{environment}}
#' @export
#' @examples      
#' evalFunctionOnList(function() rnorm(n,mean,sd),list(n=5,mean=5,sd=1))
evalFunctionOnList <-
function(x,variables=list(),...,parent_env){
    e <- new.env()
		if (!missing(parent_env))
			parent.env(e) <- parent_env
    mapply(assign, MoreArgs = list(envir=e),x=names(variables),value=variables)
    environment(x)<-e
    x(...)
}
