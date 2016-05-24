# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 10 dec 2009
# licence	: GPL

# Author(s) : Juhui Wang,Hervé Monod, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 148                 $: revision number of the last spread
#' $Author:: hmonod           $: author of the last spread
#' $Date:: 2011-06-05 23:02:3#$: date of the last spread

#-----------------------------------------------------------------------
#' Nota:
#'  A mtkProcess subclass used to perform the default analysis associated
#'  with a given sampling method (basically act as the "tell" function defined in
#'  the "sensitivity" package).  
#'  @title The mtkDefaultAnalyser class
#'  @exportClass mtkDefaultAnalyser
#=======================================================================

setClass("mtkDefaultAnalyser",
		contains=c("mtkAnalyser")
)

#' The constructor of the class \code{\linkS4class{mtkDefaultAnalyser}}
#'  @param parameters a list of parameters associated with the default Analyser method
#'  @return an object of class \code{\linkS4class{mtkDefaultAnalyser}}
#'  @examples mtkDefaultAnalyser()
#'  @export mtkDefaultAnalyser
#' @title The constructor

mtkDefaultAnalyser <- function() {
	res <- new("mtkDefaultAnalyser")


	return(res)
}

#' Performs the sensitivity analysis associated by default with the sampling method given in
#'  the \code{\linkS4class{mtkExpWorkflow}} context.
#' @title The run method 
#' @param this an object of class \code{\linkS4class{mtkDefaultAnalyser}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run

setMethod(f="run", signature=c(this="mtkDefaultAnalyser",
				context="mtkExpWorkflow"),
		definition=function(this, context){
			
			nameThis<-deparse(substitute(this)) 
			# "factors" : un objet de classe mtkExpFactors
			# "response" : un objet de classe mtkEvaluatorResult
			# Preliminary stage
			method <- context@processesVector$design@service 
			cModel<- sprintf("mtk%sAnalyser",method)
			if (!exists(cModel)) stop(paste("The class", cModel, "was not defined in this context!", sep=" "))
			analyseur<-new(cModel,name="analyze", protocol="R", site="mtk", service=method,
					parameters=context@processesVector$design@parameters, ready=TRUE, state=FALSE)
			run(analyseur,context)
			this<-analyseur
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		}
)

