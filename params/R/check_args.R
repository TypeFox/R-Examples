#' Assert none of the arguemnts of a function are null.
#'
#' Checks all the arguments in the parent function and makes sure that none of them
#' are NULL
#'
#'
#' @param ignore optionally ignore a few variables for checking [character vector].
#' @param select optionally only check a few variables of the function [character vector].
#'
#' @examples
#'
#' myfunc <- function(verbose = get_opts("verbose"), b = get_opts("b")){
#'   check_args()
#' }
#'
#' set_opts(verbose = 1)
#' ## this will throw an error, suggesting b is not defined properly
#' try(myfunc())
#'
#'
#' @export
check_args <- function(ignore, select){
	fn = sys.call(sys.parent())[1]
	env = parent.frame()
	args = ls(env)

	if(!missing(ignore))
		args = args[!args %in% ignore]

	if(!missing(select))
		args = args[args %in% select]

	miss = sapply(args, function(var){
		val = get(var, env)
		if(is.null(val))
			return(var)
		else
			return(NULL)
	})
	miss = unlist(miss) ## vars which are missing

	if ( !is.null(miss) ){
		message("Checking arguments for function: ", fn, "\n")
		message("value of following variables is null: '", paste(names(miss), collapse = ", "))
		stop("There are several options to fix this:
				 1. Use set_opts(va1 = 'value', var2 = 'value') format to define these variables.
				 2. If this function was called directly, you may simply supply these arguments to this function.
				 3. Add these parameters to configuration files, and then use load_opts('my.conf')")
	}
}
