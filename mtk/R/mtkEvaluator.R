# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 8 feb 2011
# licence	: GPL

# Author(s) : Juhui Wang, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 285                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2012-04-23 17:16:3#$: date of the last spread

#---------------------------------------------------------------------------
# 
DEBUG.mtkEvaluator=FALSE



#' The mtkEvaluator class is a sub-class of the class \code{\linkS4class{mtkProcess}} to manage
#' the model simulation task.
#' @title The mtkEvaluator class
#' @exportClass mtkEvaluator

setClass(Class="mtkEvaluator",
		contains=c("mtkProcess"),
		prototype=prototype(name="evaluate", ready=TRUE, state=FALSE),
		validity = function(object){ object@name=="evaluate"}

)

###########################################################################
## Methods definition
###########################################################################

#' The constructor
#' @param protocol a string from "http", "system", "R" respectively representing if the process is implemented remotety, locally or as R function.
#' @param site the HTTP server where the process is implemented if remotely or the package where the process is implemented if as a R function. 
#' @param service the command or R function that implements the process.
#' @param parameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the process.
#' @param parametersList a named list representing the parameters necessary to run the process. This gives another way to specify the parameters.
#' @param ready a logical to indicate if the process is ready to run.
#' @param state a logical to indicate if the simulation data are available.
#' @param result an object of a class derived from  [\code{\linkS4class{mtkEvaluatorResult}}] to hold the results produced by the process.
#' @return an object of class \code{\linkS4class{mtkEvaluator}}
#' @export mtkEvaluator

mtkEvaluator= function(protocol="R", site="mtk", service="", parameters=NULL, parametersList=NULL, ready=TRUE, state=FALSE, result=NULL) {
	
	if(!is.null(parametersList)) parameters <- make.mtkParameterList(parametersList)
	res <- new("mtkEvaluator",protocol=protocol, site=site, service=service,
			parameters=parameters,ready=ready, state=state, result=result)
	return(res) 
}


###########################################################################

#' Launches the model simulation. 
#' @param this the underlying object of class \code{\linkS4class{mtkEvaluator}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
#' @title The run method


setMethod(f="run", signature=c(this="mtkEvaluator",context="mtkExpWorkflow"), definition=function(this, context){
			
			if(this@state) return(invisible())
			if(!is.ready(this)) stop("The context of the process has not been prepared.\n")
			if(this@service=="") stop("Without results produced off-line, you must provide a service to produce the result!")
			
			nameThis=deparse(substitute(this))
			
			if(this@protocol == "system"){
				simulateur<-mtkSystemEvaluator(service=this@service,mtkParameters=this@parameters)
			}
				
			#cas trivial du package mtk
			else {
			cModel<- sprintf("mtk%sEvaluator",this@service)
			if (!exists(cModel)) stop(paste("The class", cModel, "was not defined in this context!", sep=" "))
		
			simulateur<-new(cModel, protocol="R", site="mtk", service=this@service,
					parameters=this@parameters)
			
			}
			run(simulateur,context)
			this<-simulateur
		
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
			
			
		})
