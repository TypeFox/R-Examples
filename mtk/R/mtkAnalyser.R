# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 9 feb 2011
# licence	: GPL

# Author(s) : Juhui Wang, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 285                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2012-04-23 17:16:3#$: date of the last spread

#-----------------------------------------------------------------------

#=======================================================================
## mtkAnalyser class,  definition and methods
###########################################################################

###########################################################################
## mtkAnalyser is a class derived from the  class mtkProcess. It is used 
## to generate the experiment design

DEBUG.mtkAnalyser=FALSE


#' The class mtkAnalyser is a sub-class of \code{\linkS4class{mtkProcess}} used to manage the analysis task used in a sensitivity analysis session
#' @exportClass mtkAnalyser
#' @title The mtkAnalyser class

setClass(Class="mtkAnalyser",
		contains=c("mtkProcess"),
		prototype=prototype(name="analyze",ready=TRUE, state=FALSE),
		validity = function(object){ object@name=="analyze"}

)

###########################################################################
## Methods definition
###########################################################################

###########################################################################


#' The constructor of class [\code{\linkS4class{mtkAnalyser}}].

#' @param protocol a string from "http", "system", "R" respectively representing if the Analyser is implemented remotety, locally or as R function.
#' @param site the server where the Analyser is implemented if remotely or the package where the Analyser is implemented if as a R function. 
#' @param service the command or function that implements the Analyser.
#' @param parameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Analyser.
#' @param parametersList a named list representing the parameters necessary to run the process. This gives another way to specify the parameters.
#' @param ready a logical to indicate if the process is ready to run.
#' @param result an object of a class derived from  [\code{\linkS4class{mtkAnalyserResult}}] to hold the results produced by the Analyser.
#' @return an object of class \code{\linkS4class{mtkAnalyser}}
#' @export mtkAnalyser
#' @title The constructor 


mtkAnalyser= function(protocol="R", site="mtk", service="", parameters=NULL, parametersList=NULL, ready=TRUE, state=FALSE, result=NULL) {
	if(!is.null(parametersList)) parameters <- make.mtkParameterList(parametersList)
	res <- new("mtkAnalyser", protocol=protocol, site=site, service=service,
			parameters=parameters,ready=ready, state=state, result=result)
	return(res) 
}


###########################################################################


#' Launches the sensitivity analysis task.
#' @exportMethod run
#' @param this an object of class [\code{\linkS4class{mtkAnalyser}}].
#' @param context an object of class [\code{\linkS4class{mtkExpWorkflow}}]
#' @title The method run

setMethod(f="run", signature=c(this="mtkAnalyser",context="mtkExpWorkflow"), definition=function(this, context){
			if(this@state) return(invisible())
			if(!is.ready(this)) stop("The context of the process has not been prepared.\n")
			if(this@service=="") stop("Without results produced off-line, you must provide a service to produce the result!")
			nameThis=deparse(substitute(this))
			
			#cas trivial du package mtk
			cModel<- sprintf("mtk%sAnalyser",this@service)
			if (!exists(cModel)) stop(paste("The class", cModel, "was not defined in this context!", sep=" "))
			analyseur<-new(cModel,name="analyze", protocol="R", site="mtk", service=this@service,
					parameters=this@parameters)
			run(analyseur,context)
			this<-analyseur
			
		assign(nameThis, this, envir=parent.frame())
		return(invisible())

})
