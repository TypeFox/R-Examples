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
#' $Rev:: 257                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2011-12-07 12:05:4#$: date of the last spread

#---------------------------------------------------------------------------
# nota :
# mtkResult class and its main subclasses:
#  - mtkDesignerResult
#  - mtkEvaluatorResult
#  - mtkAnalyserResult
#---------------------------------------------------------------------------
# Class: mtkResult
#---------------------------------------------------------------------------

#'  A general and simple class to collect the results produced by diverse process. 
#'  @title The mtkResult class
#'  @slot information a list  to provide information about the context of the process.
#'  @exportClass mtkResult

setClass(Class="mtkResult",
		representation=representation(
				information="ANY"
		),
		prototype=prototype(information=list()),
)

#' The constructor
#' @param information a list  to provide information about the context of the process (parameters, attributes, etc.).
#' @return an object of class \code{\linkS4class{mtkResult}}
#' @export mtkResult

mtkResult <- function(information=list()) {
	res <- new("mtkResult", information=information)
	return(res)
}

###########################################################################

#' Shows a summary of the results.
#' @param object the underlying object of class \code{\linkS4class{mtkResult}}
#' @return a summary report of the results managed by the \code{\linkS4class{mtkResult}}.
#' @exportMethod  summary
#' @title The summary method
setMethod(f="summary", signature="mtkResult", definition=function(object,...){
			if(length(object@information)>0) summary(object@information,...)
			else cat("       none")
			
		}
)

#' Serialize the process.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a summary report of the structure of the process.
#' @exportMethod  serializeOn
#' @title The serializeOnt method
setMethod(f="serializeOn", signature=c("mtkResult"), definition=function(this){
			
			return(list(information=this@information))
		}
)
