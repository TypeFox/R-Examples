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
#' $Rev:: 299                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2012-08-23 09:51:3#$: date of the last spread

#-----------------------------------------------------------------------
#
#=======================================================================
# mtkExpWorkflow class,  definition and methods
#
###########################################################################
DEBUG.mtkExpWorkflow=FALSE



#' The class mtkExpWorkflow is the main class used to manage  a sensitivity analysis session.
#' A session of sensitivity analysis usually  manages three processes: experimental design, model simulation and sensitivity analysis.
#' @title The mtkExpWorkflow class
#' @slot expFactors NULL or an object of class \code{\linkS4class{mtkExpFactors}}.
#' 
#' @slot processesVector NULL or a vector of \code{\linkS4class{mtkProcess}}  containing the processes managed within the workflow.
#' 
#' @exportClass mtkExpWorkflow

setClass(Class="mtkExpWorkflow",
		
		representation=representation(
				expFactors="ANY",
				processesVector="ANY"
		),
		
		prototype=prototype(expFactors=NULL,processesVector=NULL)	

)




###########################################################################
## Methods definition
###########################################################################


#' The constructor
#' @param expFactors an object of class \code{\linkS4class{mtkExpFactors}}.
#' @param processesVector a vector of objects from \code{\linkS4class{mtkProcess}} or its sub-classes.
#' @return an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @export mtkExpWorkflow

mtkExpWorkflow = function(expFactors=NULL, processesVector=NULL, xmlFilePath=NULL) {
	res <- new("mtkExpWorkflow",expFactors=expFactors,processesVector=processesVector)
	if(is.null(xmlFilePath)) return(res)
	parsor <- mtkParsor(xmlFilePath)
#'	show(parsor)
	run(parsor, res)
	return(res)
	
}


#' Assigns a new process to the workflow 
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' 
#' @param p an object of class \code{\linkS4class{mtkProcess}} or its sub-classes.
#' 
#' @param name a string from "design", "evaluate", or "analyze" to charaterize the process.
#' 
#' @title The addProcess method
#' @exportMethod addProcess

setMethod(f="setProcess", signature=c(this="mtkExpWorkflow",p="mtkProcess", name="character"),
		definition=function(this, p, name) {
			nameThis<-deparse(substitute(this))
			nom <- names(this@processesVector)
			v<-this@processesVector[name!=nom]
			nom<-names(v)
			this@processesVector<-c(v,p)
			names(this@processesVector) <- c(nom,name)
			reevaluate(this, name)
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})


#' Adds a process to the workflow 
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' @param p an object of class \code{\linkS4class{mtkProcess}} or its sub-classes.
#' @param name a string from "design", "evaluate", or "analyze" used to characterize the process.
#' @title The addProcess method
#' @exportMethod addProcess

setMethod(f="addProcess", signature=c(this="mtkExpWorkflow",p="mtkProcess", name="character"),
		definition=function(this, p, name) {
			nameThis<-deparse(substitute(this))
			nom <- names(this@processesVector)
			this@processesVector <-c(this@processesVector,p)
			names(this@processesVector) <- c(nom,name)
			reevaluate(this, name)
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})


#' Deletes a process from  the workflow 
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' @param name  a string from "design", "evaluate", or "analyze"  indicating the  process to remove.
#' @exportMethod deleteProcess
#' @title The deleteProcess method

setMethod(f="deleteProcess", signature=c(this="mtkExpWorkflow",name="character"),
		definition=function(this, name) {
			nameThis<-deparse(substitute(this))
			
			nom <- names(this@processesVector)
			this@processesVector<-this@processesVector[name!=nom]
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

#' gets a process from  the workflow 
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' @param name  a string from "design", "evaluate", or "analyze"  indicating the  process to get.
#' @exportMethod getProcess
#' @title The getProcess method

setMethod(f="getProcess", signature=c(this="mtkExpWorkflow",name="character"),
		definition=function(this, name) {
			nom <- names(this@processesVector)
			if(name=="design") return(this@processesVector$design)
			if(name=="evaluate") return(this@processesVector$evaluate)
			return(this@processesVector$analyze)
		})

#' extracts the results produced par  the workflow 
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' @param name  a vector of string from "design", "evaluate", or "analyze"  indicating the results to get.
#' @exportMethod extractData
#' @title The extractData method

setMethod(f="extractData", signature=c(this="mtkExpWorkflow",name="character"),
		definition=function(this, name) {
			d<-NULL
			for(x in name){
				if(is.null(d)) d<-getData(getProcess(this,x))
				else d<-cbind(d,getData(getProcess(this,x)))
			}
			return(d)
		})

#' Serialize the workflow.
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return a summary report of the structure of the process.
#' @exportMethod  serializeOn
#' @title The serializeOnt method
setMethod(f="serializeOn", signature=c(this="mtkExpWorkflow"), definition=function(this){
	apply(this@processesVector, serializeOn)
})
#' Launches the workflow 
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' @param context missing.
#' @title The run method
#' @exportMethod run

setMethod(f="run", signature=c(this="mtkExpWorkflow", context="missing"), 
		definition=function(this,context){
			nameThis<-deparse(substitute(this))
			
#		processes=vector(mode="raw", length=0)
#		noms=names(this@processesVector)
#		for(x in this@processesVector){ ## les process en aval ne peuvent pas utiliser les résultats produits par les processes en amont.
#			if(is.ready(x)) run(x, this)
#			processes=c(processes,x)
#		}
#		this@processesVector=processes
#		names(this@processesVector)=noms
			
			p <- this@processesVector$design
			if(! is.null(p) && is.ready(p)){
				run(p, this)
				this@processesVector$design <- p
				cat("\n step 1: Sampled with ", p@service, " method \n" )
			}
			p <- this@processesVector$evaluate
			if(! is.null(p) && is.ready(p)){
				run(p, this)
				this@processesVector$evaluate <- p
				
				cat("\n step 2: Simulated with ", p@service, " method \n" )
			}	
			p <- this@processesVector$analyze
			if(! is.null(p) && is.ready(p)) {
				run(p, this)
				this@processesVector$analyze <- p
				
				cat("\n step 3: Analyzed with ", p@service, " method \n" )
			}
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

#' Shows a summary of the results produced by the workflow.
#' @param object the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' @title The summary method
#' @exportMethod summary
		
		setMethod(f="summary", signature="mtkExpWorkflow", definition=function(object, ...)
				{
					#sapply(object@processesVector, summary, ...)
					p <- object@processesVector
					if(! is.null(p) && !is.null(p$design)) summary(p$design, ...)
					if(! is.null(p) && !is.null(p$evaluate)) summary(p$evaluate, ... )
					if(! is.null(p) && !is.null(p$analyze)) summary(p$analyze, ...)
					cat(" ")
					return(invisible())
				}
		)

		
###########################################################################
		
#' Show a graphical plot of the results produced by the workflow.
#' @param x the underlying object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return a graphical report of the results produced by the \code{\linkS4class{mtkExpWorkflow}}.
#' @exportMethod  plot
#' @title The plot method
		setMethod(f="plot", signature="mtkExpWorkflow", definition=function(x,y,...){
					
					if(missing(y))sapply(x@processesVector, plot, ...) 
					else sapply(x@processesVector, plot, y, ...)
					return(invisible())
				}
		)
		
###########################################################################
		
#' Prints the results produced by the workflow.
#' @param x the underlying object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return a summary of the results produced by the \code{\linkS4class{mtkExpWorkflow}}.
#' @exportMethod  print
#' @title The print method
		setMethod(f="print", signature="mtkExpWorkflow", definition=function(x,...){
					sapply(x@processesVector, print,...)
					return(invisible())
				}
		)

###########################################################################
#' Reports the results produced by the workflow.
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' 
#' @title The report method
#' @exportMethod report

setMethod(f="report", signature=c(this="mtkExpWorkflow"), definition=function(this)
		{
			sapply(this@processesVector, report)
			cat(" ")
			return(invisible())
		}
)

#' After changing a processus, re-evaluate  other processes of  the workflow to know if they should re-run. 
#' @param this the underlying object of class \code{\linkS4class{mtkExpWorkflow}}.
#' @param name a string from "design", "evaluate", or "analyze" used to characterize the process from which
#' the workflow needs to re-evaluate.
#' @title The reevaluate method

setMethod(f="reevaluate", signature=c(this="mtkExpWorkflow", name="character"),
		definition=function(this, name) {
			nameThis<-deparse(substitute(this))
			if(name=="analyze") return(invisible())
			if(name=="evaluate" && !is.null(this@processesVector$analyze)) this@processesVector$analyze@state<-FALSE
			if(name=="design"){
				if(!is.null(this@processesVector$evaluate))this@processesVector$evaluate@state<-FALSE
				if(!is.null(this@processesVector$analyze))this@processesVector$analyze@state<-FALSE
			}
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

