# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 19 avril 2012
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
# mtkProcess class,  definition and methods
############################################################################
DEBUG.mtkProcess=FALSE


#' The class mtkProcess is the root class of the processes managed within the workflow
#' @slot name a string to name the nature of  the processing: design, evaluate, analyze.
#' @slot protocol the protocol used to run the process: http, system, R, etc.
#' @slot site the package or the HTTP server where the task to run is located.
#' @slot service the function to call to fulfill the task
#' @slot parameters a vector of \code{\linkS4class{mtkParameter}} containing the parameters used to call the function defined in the service.
#' @slot ready a logical to tell if the process is ready to run.
#' @slot state a logical to tell  if the results produced by the process are availables.
#' @slot result  a data holder to keep the results produced by the process 
#' @exportClass mtkProcess
#' @title The mtkProcess class

setClass(Class="mtkProcess",
		
		representation=representation(
				name="character",
				protocol="character",
				site="character",
				service="character",
				parameters="ANY", 
				ready="logical",
				state="logical",
				result="ANY"
		),
		
		prototype=prototype(parameters=NULL, ready=FALSE, state=FALSE, result=NULL)

)

###########################################################################


#' The constructor
#' @param name the processing step associated with this process. It may be "design", "evaluate", or "analyze".
#' @param protocol a string from "http", "system", "R" respectively representing if the process is implemented remotety, locally or as R function .
#' @param site the HTTP server  where the process is implemented if remotely or the package where the process is implemented if as a R function. 
#' @param service the system command, the HTTP service  or the R function to call at the time of launching the process.
#' @param parameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the process.
#' @param ready a logical to indicate if the process is ready to run.
#' @param state a logical to indicate if the results produced by the process are already available. In this case, the process does not need to run.
#' @param result an object of a class derived from  [\code{\linkS4class{mtkResult}}] to hold the results produced by the process.
#' @return an object of class \code{\linkS4class{mtkProcess}}
#' @export mtkProcess
mtkProcess = function(name, protocol="R", site="mtk", service="", parameters=NULL, ready=FALSE, state=FALSE, result=NULL) {
	res <- new("mtkProcess",name=name, protocol=protocol, site=site, service=service,
			parameters=parameters,ready=ready,state=state, result=NULL)
	return(res)
	
}

###########################################################################
#' Gives a name to the process
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @param name a string indicating the processing step associated with this process. It may be "design", "evaluate", or "analyze".
#' @return invisble()
#' @exportMethod setName
#' @title The setName method

setMethod(f="setName", signature=c(this="mtkProcess", name="character"),
		definition=function(this, name ) {
			nameThis <- deparse(substitute(this))
			this@name <-  name
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

###########################################################################
#' get the name of the process
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a string indicating the processing step associated with this process.
#' @exportMethod getName
#' @title The getName method

setMethod(f="getName", signature=c(this="mtkProcess"),
		definition=function(this) {
			return(this@name)
		})

###########################################################################


#' Assigns a new vector of parameters to  the process.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @param f  a vector of  \code{\linkS4class{mtkParameter}}.
#' @return invisble()
#' @exportMethod setParameters
#' @title The setParameters method

setMethod(f="setParameters", signature=c(this="mtkProcess",f="vector"),
		definition=function(this,f) {
			nameThis <- deparse(substitute(this))
			this@parameters <- f
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

###########################################################################

#' Returns the vector of parameters as a named list.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a named list in which each element corresponds to a parameter. The vector of parameters are represented as  a named list such as 
#' (name of parameter1 = value of parameter1, name of parameter2 = value of parameter 2, ...).
#' @exportMethod getParameters
#' @title The getParameters method

setMethod(f="getParameters",
		signature=c(this="mtkProcess"),
		definition=function(this) {
			if(is.null(this@parameters)) return(NULL)
			n<-lapply(this@parameters,getName)
			v<-lapply(this@parameters,getValue)
			names(v)=n
			return(v)
		}
)

###########################################################################

#' Returns the result.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return an object of class \code{\linkS4class{mtkResult}}
#' @exportMethod getResult
#' @title The getResult method

setMethod(f="getResult",
		signature=c(this="mtkProcess"),
		definition=function(this) {
			return(this@result)
		}
)

###########################################################################

#' Returns the result as a data frame.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a data frame corresponding the attribut "main" of the class \code{\linkS4class{mtkResult}}
#' @exportMethod getData
#' @title The getData method

setMethod(f="getData",
		signature=c(this="mtkProcess"),
		definition=function(this) {
			return(this@result@main)
		}
)
###########################################################################


#' Tests if the process is ready to start.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return TRUE or FALSE.
#' @exportMethod is.ready
#' @title The is.ready method

setMethod(f="is.ready", signature=c("mtkProcess"), 
		definition=function(this) {
			return(this@ready)
		})

###########################################################################
#' Tests if the results produced by the process are available.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return TRUE or FALSE.
#' @exportMethod is.finished
#' @title The is.finished method

setMethod(f="is.finished", signature=c("mtkProcess"), 
		definition=function(this) {
			return(this@state)
		})

###########################################################################


#' Makes the process ready to run.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @param switch  a logical (TRUE or FALSE). 
#' @return invisble()
#' @exportMethod setReady
#' @title The setReady method

setMethod(f="setReady", signature=c(this="mtkProcess",switch="logical"),
		definition=function(this, switch) {
			nameThis <- deparse(substitute(this))  
			this@ready <- switch
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
###########################################################################


#' Tells the process that the results that it should produce are available.
#' Such results might be produced off-line without running the process.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @param state  a logical (TRUE or FALSE). 
#' @return invisble()
#' @exportMethod setState
#' @title The setState method

setMethod(f="setState", signature=c(this="mtkProcess",state="logical"),
		definition=function(this, state) {
			nameThis <- deparse(substitute(this))  
			this@state <- state
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

###########################################################################



#' Serialize the process.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a summary report of the structure of the process.
#' @exportMethod  serializeOn
#' @title The serializeOnt method
setMethod(f="serializeOn", signature=c("mtkProcess"), definition=function(this){
			parametresListe<-NULL
			if(!is.null(this@parameters))
				parametresListe<-getParameters(this)## parametresList<- serialisation de la list des parametres
			
			resultat<-NULL
			if(!is.null(this@result)) resultat<-serializeOn(this@result)
				
			return(list(attributList=
							list(name=this@name,
							protocol=this@protocol,
							site=this@site,
							service=this@service,
							parameters=parametresListe, 
							ready=this@ready,
							state=this@state,
							result=resultat),
							classname="mtkProcess"))
		}
)

###########################################################################

#' Shows a summary of the results.
#' @param object the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a summary report of the results produced by the process.
#' @exportMethod  summary
#' @title The summary method
setMethod(f="summary", signature=c("mtkProcess"), definition=function(object,...){
		
	
		cat("---------------------------------------------------------------------------\n")
		
		cat(paste("	 HERE IS A SUMMARY OF THE PROCESS \"", object@name, "\": \n", sep=""))
		cat("---------------------------------------------------------------------------\n\n\n")
		cat("TYPE : ", object@name, "\n")
		cat("METHOD : ", object@service, "\n\n")
		if(object@state) summary(object@result,...)
		cat("\n\n")
	}
		
)

###########################################################################

#' Plots the results produced by the process.
#' @param x the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a graphical report of the results produced by the process.
#' @exportMethod  plot
#' @title The plot method
setMethod(f="plot", signature=c("mtkProcess"), definition=function(x,y,...){
			if(x@state) { 
				dev.new()
				if(!missing(y)) plot(x@result,y,...)
				else plot(x@result,...)
			}
			
		}
)

###########################################################################

#' Prints the results produced by the process.
#' @param x the underlying object of class \code{\linkS4class{mtkProcess}}.
#' @return a summary of the results produced by the \code{\linkS4class{mtkProcess}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="print", signature=c("mtkProcess"), definition=function(x,...){
			cat("---------------------------------------------------------------------------\n")
			cat(paste("	 HERE IS INFORMATION ABOUT THE PROCESS \"", x@name, "\": \n", sep=""))
			cat("---------------------------------------------------------------------------\n\n\n")
			cat("TYPE : ", x@name, "\n")
			cat("METHOD : ", x@service, "\n\n")
			if(x@state) print(x@result,...)
			invisible()
		}
)

###########################################################################
#' Reports the results.
#' @param this the underlying object of class \code{\linkS4class{mtkProcess}}
#' @return a summary report of the results produced by the process.
#' @exportMethod  report
#' @title The report method
setMethod(f="report", signature=c(this="mtkProcess"), definition=function(this){
			cat("-------------------------------------------------------\n")
			cat(paste("--------- THE RESULTS OF THE \"", this@name, "\" ARE ---------------\n\n\n", sep=""))
			if(this@state) show(this@result)
			cat("\n\n")
		}
)
