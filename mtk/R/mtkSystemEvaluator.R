# Mexico Toolkit
# Author(s) : J. WANG, INRA-MIA-Jouy, 78352, Jouy en Josas, France 
# Repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/


#' A sub-class of the class \code{\linkS4class{mtkEvaluator}} used to perform model simulation
#' with a model implemented with a external application.
#' @title The mtkSystemEvaluator class
#' @exportClass mtSystemEvaluator

setClass("mtkSystemEvaluator",
		contains=c("mtkEvaluator")
	)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the evaluator.
#' @param listParameters a named list defining the parameters necessary to run the evaluator. It gives non object-oriented way to define the parameters.
#' @return an object of class \code{\linkS4class{mtkSystemEvaluator}}
#' @examples mtkSystemEvaluator()
#' @export mtkSystemEvaluator
#' @title The constructor

mtkSystemEvaluator<- function(service="", mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
		if(!is.null(listParameters))
			p <- make.mtkParameterList(listParameters)
					
	res <- new("mtkSystemEvaluator",protocol="system", service=service, parameters=p)
		return(res)
		}


#' Performs the simulation  with a model implemented as an external application. 
#' @title The run method
#' @param this an object of class \code{\linkS4class{mtkSystemEvaluator}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
	
setMethod(f="run", signature=c(this="mtkSystemEvaluator",
			context="mtkExpWorkflow"),
		definition=function(this, context){
			if(this@state) return(invisible())
			nameThis<-deparse(substitute(this))
			
			X <- context@processesVector$design@result@main 
			parameters<-getParameters(this)
			
	
			output<-system(this@service, intern=FALSE,wait=TRUE )
			if(output==0) show("Tout se passe bien avec l'appel system")
	
			this@result <- mtkSystemEvaluatorResult(main=data.frame())
			this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
	
