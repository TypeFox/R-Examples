# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 10 dec 2009
# licence	: GPL

# Author(s) : Hervé Monod, Juhui Wang, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 269                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2012-01-17 10:41:3#$: date of the last spread
#
#---------------------------------------------------------------------------
# Class: mtkEvaluatorResult
#---------------------------------------------------------------------------

#'  A class to collect simulation results 
#'  @title The mtkEvaluatorResult class
#'  @slot main a data-frame containing the simulation results
#'  @exportClass mtkEvaluatorResult


setClass(Class="mtkEvaluatorResult",
		representation=representation(
				main="data.frame"
		),
		contains=c("mtkResult"),
		prototype=prototype(main=data.frame()),
)

#' The constructor
#' @param main a data-frame to hold the results of simulation.
#' @param information a named list to give supplementary information about the main data andthe sampling process.
#' @return an object of class \code{\linkS4class{mtkEvaluatorResult}}
#' @export mtkEvaluatorResult

mtkEvaluatorResult <- function(main=data.frame(),
		information=list()
){
	res <- new("mtkEvaluatorResult", main=main, information=information)
	return(res)
}

##########################################################################

#' Shows a summary of the results.
#' @param object the underlying object of class \code{\linkS4class{mtkEvaluatorResult}}
#' @return a summary report of the results managed by the \code{\linkS4class{mtkEvaluatorResult}}.
#' @exportMethod  summary
#' @title The summary method
setMethod(f="summary", signature="mtkEvaluatorResult", definition=function(object,...){
			cat("ABOUT USED PARAMETERS : \n")
			if(length(object@information)>0) show(summary(object@information,...))
			else cat("      none")
			cat("\n")
			cat("ABOUT SIMULATION RESULTS :\n\n")
			cat("      number of simulations :",nrow(object@main),"\n\n")
			 show(summary(object@main,...))
			
		}
)

###########################################################################

#' Show a graphical plot of the results.
#' @param x the underlying object of class \code{\linkS4class{mtkEvaluatorResult}}.
#' @return a graphical report of the results managed by the \code{\linkS4class{mtkEvaluatorResult}}.
#' @exportMethod  plot
#' @title The plot method
setMethod(f="plot", signature="mtkEvaluatorResult", definition=function(x, y, ...){
		if(missing(y))plot(x@main,...)
					else plot(x@main, y, ...)
			
		}
)

###########################################################################

#' Prints the managed data.
#' @param x the underlying object of class \code{\linkS4class{mtkEvaluatorResult}}.
#' @return a summary of the results managed by  the class \code{\linkS4class{mtkEvaluatorResult}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="print", signature="mtkEvaluatorResult", definition=function(x,...){
			
				cat("USED PARAMETERS : \n")
				if(length(x@information)>0) print(x@information,...)
				else cat("       none")
				cat("\n")
				cat("SIMULATION RESULTS :\n\n")
			
				print(x@main,...)
			
			
		}
)