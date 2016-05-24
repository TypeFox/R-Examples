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
#' $Rev:: 142                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2011-06-01 17:31:2#$: date of the last spread

#-----------------------------------------------------------------------

#=======================================================================
# Class: mtkFastDesigner
#---------------------------------------------------------------------------
#'  A mtkProcess subclass to produce the experiment design  with the Fast method
#'  @title The mtkFastDesigner class
#'  @exportClass mtkFastDesigner

setClass(Class="mtkFastDesigner",
		contains=c("mtkDesigner")
)

#' The constructor.
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Designer.
#' @param listParameters a named list defining the parameters necessary to run the Designer.
#' @return an object of class \code{\linkS4class{mtkFastDesigner}}
#' @examples mtkFastDesigner()
#' @export mtkFastDesigner
#' @title The constructor

mtkFastDesigner <- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
	if(!is.null(listParameters))
		p <- make.mtkParameterList(listParameters)
	res <- new("mtkFastDesigner",service="Fast", parameters=p)
	return(res)
}

#' Builds the experiment design with  the Fast method
#' @param this the underlying object of class \code{\linkS4class{mtkFastDesigner}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @note based on the fast99 function of the "sensitivity" package
#' @exportMethod run
#' @title  The run method 
setMethod(f="run", signature=c(this="mtkFastDesigner",
				context="mtkExpWorkflow"),
		definition=function(this, context){
			if(this@state) return(invisible())
			nameThis = deparse(substitute(this))
			
			context=context@expFactors ## Line added by Juhui
			
			factorNames <- getNames(context)
			distribNames <- getDistributionNames(context)
			distribParameters <- getDistributionParameters(context)
			processParameters <- getParameters(this)
			
			##!!
			##!! Pre-processing the input data, the processing of the method to implement follows:
			##!!
			
			# Calculations: call to the fast99 function of the sensitivity library
			resultat <- eval( do.call("fast99",
			c(list(factors=factorNames,
                               q=paste("q",distribNames,sep=""),
                               q.arg=distribParameters),
                          processParameters))
			)
			# Output
			information <- resultat[c("model", "M", "s", "omega", "call")]
			information$SamplingMethod <- "fast99"
			design <- as.data.frame(resultat$X)
			
			##!!
			##!! post-processing the output of the method:
			##!!
			
			this@result <- mtkDesignerResult(main=design, information=information)
			this@state=TRUE
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		}
)
#' A sub-class of the class \code{\linkS4class{mtkDesignerResult}} used to hold the results of the sensitivity  analysis
#' @title The mtkFastDesignerResult class
#' @exportClass mtkFastDesignerResult

setClass("mtkFastDesignerResult",
		
		contains=c("mtkDesignerResult")
)

#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Designer.
#'  @param information a named list to provide supplementary information about the analysis process and its results.

#' @return an object of class \code{\linkS4class{mtkFastDesignerResult}}
#' @examples mtkmtkFastDesignerResult()
#' @export mtkmtkFastDesignerResult
#' @title The constructor

mtkFastDesignerResult <- function(main, information=NULL) {
	res <- new("mtkFastDesignerResult", main=main, information=information)
	return(res)
}