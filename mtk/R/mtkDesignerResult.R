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
# Class: mtkDesignerResult
#---------------------------------------------------------------------------

#'  A class to hold the experiment design produced by a Designer.
#'  @title The mtkDesignerResult class
#'  @slot main a data-frame containing the produced experiment design.
#'  @exportClass mtkDesignerResult

setClass(Class="mtkDesignerResult",
		representation=representation(
				main="data.frame"
		),
		contains=c("mtkResult"),
		prototype(main=data.frame())
)

#' The constructor
#'  @param main a data-frame to hold the main results of the experiment design.
#'  @param information a named list to provide supplementary information about the sampling process and its results.
#'  @return an object of class \code{\linkS4class{mtkDesignerResult}}
#'  @export mtkResult
#'  @title The constructor
mtkDesignerResult <- function(main=data.frame(),
		information=list())
{
	res <- new("mtkDesignerResult", main=main, information=information)
	return(res)
}


###########################################################################

#' Shows a summary of the results.
#' @param object the underlying object of class \code{\linkS4class{mtkDesignerResult}}
#' @return a summary report of the results managed by the \code{\linkS4class{mtkDesignerResult}}.
#' @exportMethod  summary
#' @title The summary method
setMethod(f="summary", signature="mtkDesignerResult", definition=function(object,...){
			
			cat("ABOUT THE PARAMETERS USED : \n")
			if(length(object@information)>0) show(summary(object@information, ...))
			else cat("       none")
			cat("\n")
			cat("ABOUT THE GENERATED EXPERIMENT DESIGN :\n\n")
			cat("       number of simulations :",nrow(object@main),"\n\n")
			show(summary(object@main,...))
		}
)

###########################################################################

#' Show a graphical plot of the results.
#' @param x the underlying object of class \code{\linkS4class{mtkDesignerResult}}.
#' @return a graphical report of the results managed by the \code{\linkS4class{mtkDesignerResult}}.
#' @exportMethod  plot
#' @title The plot method
setMethod(f="plot", signature="mtkDesignerResult", definition=function(x,y, ...){
			if(!missing(y))plot(x@main, y, ...)
			else plot(x@main, ...)
			
		}
)

###########################################################################

#' Prints the managed data.
#' @param x the underlying object of class \code{\linkS4class{mtkDesignerResult}}.
#' @return a summary of the results managed by the \code{\linkS4class{mtkDesignerResult}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="print", signature="mtkDesignerResult", definition=function(x,...){
			
				cat("USED PARAMETERS : \n")
				if(length(x@information)>0) print(x@information,...)
				else cat("       none")
				cat("\n")
				cat("EXPERIMENTS DESIGN :\n\n")
			
				print(x@main,...)
			
			
		}
)