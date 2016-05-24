#' Interpreting Time Series, Autocorrelated Data Using GAMMs (itsadug)
#'
#' Itsadug provides a set of functions that facilitate the evaluation, 
#' interpretation, and visualization of GAMM models that are implemented in 
#' the package \code{\link[mgcv]{mgcv}}. 
#' 
#' @section Tutorials:
#' \itemize{
#' \item \code{vignette("inspect", package="itsadug")} - summarizes different 
#' functions for visualizing the model.
#' \item \code{vignette("test", package="itsadug")} - summarizes 
#' different functions for significance testing.
#' \item \code{vignette("acf", package="itsadug")} - summarizes how to check 
#' and account for autocorrelation in the residuals.
#' }
#' Also available online: \url{www.jacolienvanrij.com/itsadug}.
#'
#' @section Interpretation and visualization:
#' Main functions that are provided in \code{itsadug} for interpretation and 
#' visualization of GAMM models:
#' \itemize{
#'   \item \code{\link{pvisgam}} plots partial interaction surfaces; it also 
#' allows for visualizing 3-way or higher interactions.
#'   \item \code{\link{fvisgam}} plots summed interaction surfaces, with the 
#' possibility to exclude random effects.
#'   \item \code{\link{plot_smooth}} plots 1D model estimates, and has the 
#' possibility to exclude random effects. 
#'   \item \code{\link{plot_parametric}} plot group estimates.
#'   \item \code{\link{inspect_random}} plots and optionally averages random 
#' smooths
#'   \item \code{\link{plot_data}} plots the data
#'   \item \code{\link{plot_topo}} plots EEG topographies
#' }
#' 
#' @section Testing for significance:
#' \itemize{
#' \item \code{\link{compareML}} Performs Chisquare test on two models
#' \item \code{\link{plot_diff}} Calculates and visualizes the difference 
#' between two conditions within a model
#' \item \code{\link{plot_diff2}} Calculates and visualizes the 2 dimensional 
#' difference between two conditions within a model
#' }
#' 
#' @section Evaluation of the model:
#' \itemize{
#' \item \code{\link{check_resid}} plots four different plots to inspect the 
#' distribution of and structure in the residuals
#' \item \code{\link{plot_modelfit}} plots an overlay of the data and the 
#' modelfit for randomly selected trials
#' \item \code{\link{diagnostics}} produces plots of the distributions of 
#' residuals and predictors in the model
#' }
#' 
#' @section Checking and handling autocorrelation:
#' \itemize{
#' \item \code{\link{acf_resid}} different ways to inspect autocorrelation in 
#' the residuals
#' \item \code{\link{start_event}} creates an AR.start column 
#' \item \code{\link{resid_gam}} returns residuals corrected for the AR1 model
#' }
#'
#' @section Predictions:
#' Further, there are some wrappers around the \code{\link[mgcv]{predict.gam}}
#' function to facilitate the extraction of model predictions. These can be 
#' used for customized plots. See for an example in the vignette 
#' "plotfunctions" 
#' (\code{vignette("plotfunctions", package="itsadug")}).
#' \itemize{
#'   \item \code{\link{get_predictions}} for getting the estimates for given 
#' settings of some or all of the model predictors;
#'   \item \code{\link{get_difference}} for extracting the difference between 
#' two conditions or two smooths or two surfaces.
#'   \item \code{\link{get_modelterm}} for extracting the smooth term (
#' partial) estimates.
#'   \item \code{\link{inspect_random}} and \code{\link{get_random}} for 
#' extracting random effects only.
#' }
#'
#' @section Notes:
#' \itemize{
#' \item Use \code{\link{infoMessages}(FALSE)} to suppress all 
#' information messages for the current session. 
#' This may be helpful when creating knitr or
#' R markdown reports.
#' \item The vignettes are available via \code{browseVignettes()}. 
#' When working on a server via the command line, 
#' using \code{ssh -X} instead of \code{ssh} may make the 
#' HTML files available.
#' \item A list of all available functions is provided in 
#' \code{help(package="itsadug")}.
#' }
#'
#' @author
#' Jacolien van Rij, Martijn Wieling, R.Harald Baayen, Hedderik van Rijn
#'
#' Maintainer: Jacolien van Rij (\email{vanrij.jacolien@gmail.com})
#'
#' University of Groningen, The Netherlands & University of Tuebingen, Germany
#' @docType package
#' @name itsadug
NULL





#' Turn on or off information messages.
#' 
#' @export
#' @param input Input variable indicating to print info messages 
#' ("on", or 1, or TRUE) or not ("off", 0, or FALSE).
#' @examples
#' # To turn on the info messages (all the same):
#' infoMessages("on")
#' infoMessages(1)
#' infoMessages(TRUE)
#' # To turn off the info messages (all the same):
#' infoMessages("off")
#' infoMessages(0)
#' infoMessages(FALSE)
#' # checking output:
#' (out <- infoMessages(FALSE))
#' @family Functions for package use
infoMessages <- function(input){
	if(is.logical(input)){
		options(itsadug_print=input)
	}else if(is.numeric(input)){
		options(itsadug_print=ifelse(input<=0, FALSE, TRUE))
	}else if(is.character(input)){
		options(itsadug_print=ifelse(input=="off", FALSE, 
			ifelse(input=="on",TRUE, getOption('itsadug_print'))))
	}else{
		stop(sprintf("Cannot interpret input value %s. Try to use logical values TRUE or FALSE.", input))
	}
	invisible(list(value=getOption('itsadug_print'), 
		effect=ifelse(getOption('itsadug_print')==TRUE, "messages printed", "no messages printed")))
}





.onAttach <- function(...) {
  if(is.null(getOption('itsadug_print'))){
  	options(itsadug_print=TRUE)
  }
  if(getOption('itsadug_print')){
  	packageStartupMessage('Loaded package itsadug 2.0 (see \'help("itsadug")\' ).')
  }
}





#' Information on how to cite this package
#' 
#' @export
#' @import utils
#' @param input Optional parameter. Normally (NULL) the citation info is 
#' printed. If value "version" then only the version is printed.
#' @examples
#' info()
#' info("version")
#' citation(package="itsadug")
#' # To get info about R version:
#' R.version.string
#' @seealso
#' \code{\link[utils]{citation}}, \code{\link[base]{R.version}},
#' \code{\link[utils]{sessionInfo}}
#' @family Functions for package use
info <- function(input=NULL){
	if(is.null(input)){
		citation(package="itsadug")
	}else if(input=="version"){
		cat(sprintf("Package itsadug, version %s\n", 
			packageVersion("itsadug")))
	}else{
		help(package="itsadug")
	}
}





