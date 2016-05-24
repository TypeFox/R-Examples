#' avisSetup
#'
#' Sets up settings that apply to the behabiour of the package
#' Allow users to turn off the information messages of the functions.
#'  
#' @usage avisSetup (...)
#' @param ... Package settings parameters. Available params: verbose = TRUE/FALSE
#' @export
#' @examples \dontrun{
#' avisSetup(verbose=FALSE)
#' }
#'

avisSetup <- function(...){
	l <- list(...)

	if(is.element('verbose', names(l))){
		.setAvisVerbosity(l[['verbose']])
	}
}

.setAvisVerbosity <- function(v){
	if("logical" != typeof(v)){
		stop("Verbosity must be of 'logical' type")
	}
	.avisCacheSet(".ravis_verbose", v)
}

.isAvisVerbose <- function(){
	.avisCacheReturnOrSetup(".ravis_verbose", function(){ 
		default_verbosity = FALSE
		default_verbosity
	})
}

# Shows message depending on library configuration
.avisVerboseMessage <- function(mens){
	verb = .isAvisVerbose()

	if(verb){
		message(mens)
	}
}