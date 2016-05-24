



#' funr: providing a simple command-line interface to R functions
#'
#' @description
#' Wraps Rscript in a easy to use manner, exposing all R functions from the terminal.
#' The \href{https://github.com/sahilseth/funr}{github page} provides more details with examples,
#' highlights and caveats.
#'
#' @aliases funr rfun cli
#'
#' @param args Should always be: \code{commandArgs(trailingOnly = TRUE)}, when used
#' inside a script. \href{https://github.com/sahilseth/funr/blob/master/inst/scripts/funr}{Example}
#' @param help_text A simple text to be displayed describing options and usage of this interface.
#' Supplying this, replaces the default text.
#' @param script_name Name of the script. This is used in the the help text. [funr]
#'
#'
#' @source https://github.com/sahilseth/funr
#'
#' @export
#' @importFrom utils head help str
#'
#' @examples
#' ## show funr help
#' ## terminal version: funr -h
#' funr()
#'
#'
#' ## show rnorm help
#' ## terminal version: funr -h rnorm
#' render_funr(funr(args=c("-h", "rnorm")))
#'
#' ## Generate a few random numbers
#' ## terminal version: funr rnorm n=10
#' render_funr(funr(args=c("rnorm", "n=10")))
#'
funr <- function(args,
								 help_text,
								 script_name = "funr"){

	##        show help if there are no arguments
	#if(missing(help_text))

	if(missing(args)){
		message(generic_help(help_text = help_text, script_name = script_name))
		return()
	}

	if(length(args) == 0){
		message(generic_help(help_text = help_text, script_name = script_name))
		return()
	}

	# Arguments which start with - are for this script
	rm = grep("^-", args)
	script_args = args[rm]

	verbose = FALSE
	if("-v" %in% script_args)
		verbose = TRUE

	if(verbose){
		message("args:"); message(args)
		message("script_args:"); message(script_args)
	}

	# remove these from subsequent processing
	if(length(rm) > 0)
		args = args[-rm]

	if(length(args) == 0){
		message(generic_help(help_text = help_text))
		return()
	}

	#      Get name of the function
	func = args[1]
	#    all arguments to that function
	args = args[-1]

	if(verbose){
		message("\nusing func:");message(func)
		message("with final args:");message(args)
	}

	args_w_missing_eq = grep("=", args, value = TRUE, invert = TRUE)
	if(length(args_w_missing_eq) > 0){
		message("> All args must use the format variable=value!\n",
		 "> Please check the followings args:",
		paste(args_w_missing_eq, collapse = "\n"),
		"\n> A valid example would be say, 'rnorm n=100'")
		generic_help()
		invisible()
	}

	#           Load the required package
	if(grepl("::", func)){
		pkg <- gsub("(.?)::.*", "\\1", func)
		cat("loading pkg:", pkg, "\n");
		library(pkg, character.only = TRUE)
		func = as.character(gsub(".*::(.*)", "\\1", func))
	}

	if( is.na(func) ){
		generic_help()
		invisible()
	}

	fn = try(get(func), silent = TRUE)
	if(class(try(fn)) == "try-error")
		stop("\nI could not find a function called '", func, "' name, please check.")

	if(is.function(fn) & "-h" %in% script_args){
		out = withVisible(help(func))
		class(out) = c("funr", "list")
		return(out)

	}else{

		params <- parse_params(func = func, paramPairs = args, verbose = verbose)
		if(verbose){
			cat("\nStarting", func, "with params\n",
					paste(names(params), unlist(params),sep=": ",
								collapse="\n"),"\n")
			#message(args)
			if(verbose) message(str(params))
		}

		if(length(args) == 0)
			message("\ntry:      ", script_name, " -h ", func, "     to get more details on this function.")

		out = try(withVisible(do.call(func, args = params)))

		if(class(out) == "try-error"){
			stop("")
		}

		class(out) = c("funr", "list")
		return(out)
	}

}
