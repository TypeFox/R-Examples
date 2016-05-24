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
#' $Rev:: 266                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2012-01-10 13:59:4#$: date of the last spread

#-----------------------------------------------------------------------

#=======================================================================
#---------------------------------------------------------------------------
# Class: mtkAnalyserResult
#---------------------------------------------------------------------------

#'  A class to hold the results produced by the Analyser
#'  @slot main a data-frame containing the analysis results.
#'  @exportClass mtkAnalyserResult
#'  @title The mtkAnalyserResult class

setClass(Class="mtkAnalyserResult",
		representation=representation(
				main="ANY"
		),
		contains=c("mtkResult"),
)

#' The constructor of class \code{\linkS4class{mtkAnalyserResult}}
#'  @param main a data-frame to hold the main results produced by the Analyser.
#'  @param information a named list to provide supplementary information about the analysis process and its results.

#'  @return an object of class \code{\linkS4class{mtkAnalyserResult}}
#'  @examples mtkAnalyserResult()
#'  @export mtkAnalyserResult
#'  @title The constructor
mtkAnalyserResult <- function(main=data.frame(), information=list()) {
	res <- new("mtkAnalyserResult", main=main, information=information)
	return(res)
}

###########################################################################

#' Shows a summary of the results.
#' @param object the underlying object of class \code{\linkS4class{mtkAnalyserResult}}
#' @return a summary report of the results managed by the \code{\linkS4class{mtkAnalyserResult}}.
#' @exportMethod  summary
#' @title The summary method
setMethod(f="summary", signature="mtkAnalyserResult", definition=function(object,...){
			
			if(class(object@main)=='list'){
				print(object,...)
				return(invisible())
			}
			cat("ABOUT USED PARAMETERS : \n")
			if(length(object@information)>0) print(summary(object@information,...))
			else cat("       none")
			cat("\n")
			cat("ABOUT ANALYSIS RESULTS :\n\n")
			print(summary(object@main,...))
			
		}
)

###########################################################################

#' Show a graphical plot of the results.
#' @param x the underlying object of class \code{\linkS4class{mtkAnalyserResult}}.
#' @return a graphical report of the results managed by the \code{\linkS4class{mtkAnalyserResult}}.
#' @exportMethod  plot
#' @title The plot method
setMethod(f="plot", signature="mtkAnalyserResult", definition=function(x, y, ...){
	
			if(missing(y))plot(x@main,...)
				else plot(x@main, y, ...)
			
		}
)

###########################################################################

#' Prints the managed data.
#' @param x the underlying object of class \code{\linkS4class{mtkAnalyserResult}}.
#' @return a summary of the results managed by  the class \code{\linkS4class{mtkAnalyserResult}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="print", signature="mtkAnalyserResult", definition=function(x,...){
			
		
			cat("USED PARAMETERS : \n")
			if(length(x@information)>0) print(x@information)
			else cat("       none")
			cat("\n")
			cat("ANALYSIS RESULTS :\n\n")
		
			print(x@main,...)
			
		}
)