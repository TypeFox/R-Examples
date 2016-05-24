# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 10 dec 2009
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
##########################################################################
## mtkDesigner class,  definition and methods
###########################################################################

DEBUG.mtkDesigner=FALSE

#' The mtkDesigner class is a sub-class of the class \code{\linkS4class{mtkProcess}} used to manage 
#' the sampling task.
#' @title The mtkDesigner class
#' @exportClass mtkDesigner

setClass(Class="mtkDesigner",
		contains=c("mtkProcess"),
		prototype=prototype(name="design", ready=TRUE, state=FALSE),
		validity = function(object){ object@name=="design"}
		
)

###########################################################################


#' The constructor
#' @param protocol a string from "http", "system", "R" respectively representing if the process is implemented remotety, locally or as R function.
#' @param site the HTTP server where the process is implemented if remotely or the package where the process is implemented if as  R function. 
#' @param service the command or R function that implements the process.
#' @param parameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the process.
#' @param parametersList a named list representing the parameters necessary to run the process. This gives another way to specify the parameters.
#' @param ready a logical to indicate if the process is ready to run.
#' @param state a logical to indicate if the sampling data are produced off-line.
#' @param result an object of  class  [\code{\linkS4class{mtkDesignerResult}}] to hold the results produced by the process.
#' @return an object of class \code{\linkS4class{mtkDesigner}}
#' @export mtkDesigner

mtkDesigner= function(protocol="R", site="mtk", service="", parameters=NULL, parametersList=NULL, ready=TRUE, state=FALSE, result=NULL) {
	if(!is.null(parametersList)) parameters <- make.mtkParameterList(parametersList)
	
	res <- new("mtkDesigner", protocol=protocol, site=site, service=service,
			parameters=parameters,ready=ready, state=state, result=result)
	return(res) 
}


###########################################################################


#' Generates the design by sampling the factors. This is a virtual method and should be implemented by the derived classes.
#' @param this the underlying object of class \code{\linkS4class{mtkDesigner}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
#' @title The run method


setMethod(f="run", signature=c(this="mtkDesigner",context="mtkExpWorkflow"), definition=function(this, context){
			
			if(this@state) return(invisible())
			if(!is.ready(this)) stop("The context of the process has not been prepared.\n")
			if(this@service=="") stop("Without results produced off-line, you must provide a service to produce the result!")
			nameThis=deparse(substitute(this))
			
				#cas trivial du package mtk
				cModel<- sprintf("mtk%sDesigner",this@service)
				if (!exists(cModel)) stop(paste("The class", cModel, "was not defined in this context!", sep=" "))
				sampleur<-new(cModel,name="design", protocol="R", site="mtk", service=this@service,
						parameters=this@parameters)
				run(sampleur,context)
				this<-sampleur
				
			if(DEBUG.mtkDesigner) cat("**** Le contenu du Designer: **** \n")
			if(DEBUG.mtkDesigner) show(this)
			if(DEBUG.mtkDesigner) cat("**** Fin du contenu du Designer: *****\n\n")
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
			
			
		})

