# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 10 Feb 2011
# licence	: GPL

# Author(s) : Juhui Wang, Hervé Monod, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 176                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2011-06-21 11:11:3#$: date of the last spread

#-----------------------------------------------------------------------
#'  A sub-class of the class \code{\linkS4class{mtkAnalyser}} used to perform the  analysis 
#'  with the "fast99" method defined in the "sensitivity" package. 
#'  @title The mtkFastAnalyser class
#'  @exportClass mtkFastAnalyser
#=======================================================================
setClass("mtkFastAnalyser",
		contains=c("mtkAnalyser")
)

#' The constructor
#' @param mtkParameters a vector of [\code{\linkS4class{mtkParameter}}] representing the parameters necessary to run the Analyser.
#' @param listParameters a named list defining the parameters necessary to run the Designer. It gives non object oriented way to specify the parameters.
#'  @return an object of class \code{\linkS4class{mtkFast99Analyser}}
#'  @examples mtkFastAnalyser()
#'  @export mtkFastAnalyser
#' @title The constructor

mtkFastAnalyser <- function(mtkParameters=NULL, listParameters=NULL) {
	p<-mtkParameters
	if(!is.null(listParameters))
		p <- make.mtkParameterList(listParameters)
	res <- new("mtkFastAnalyser",service="Fast", parameters=p)
	return(res)
}

#' Performs the sensitivity analysis with the "fast99" implemented in the "sensitivity" package.
#' @title The run method 
#' @param this an object of class \code{\linkS4class{mtkFastAnalyser}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run

setMethod(f="run", signature=c(this="mtkFastAnalyser",
				context="mtkExpWorkflow"),
		definition=function(this, context){
			if(this@state) return(invisible())
			nameThis<-deparse(substitute(this)) 
			
			X <- context@processesVector$design@result@main 
			parameters<-context@processesVector$design@result@information
			Y <- context@processesVector$evaluate@result@main
			
			##!!
			##!! Pre-processing the input data, the processing of the method to implement follows:
			##!!
			
			Obj <- c(parameters,list(X=X))
			class(Obj) <- Obj$SamplingMethod
			
			
			# Management of the multivariate case (VERY BASIC at the moment)
			Y <- Y[,1]
			
			
			# Main stage:
			# here, just a call to the "tell" function of "sensitivity"
			analysisOutput <- tell(Obj, Y)
			
			##!!
			##!! post-processing the output of the method:
			##!!
			this@result <- mtkAnalyserResult(main=analysisOutput)
			this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		}
)

#' A sub-class of the class \code{\linkS4class{mtkAnalyserResult}} used to hold the results of the sensitivity  analysis
#' @title The mtkFastAnalyserResult class
#' @exportClass mtkFastAnalyserResult

setClass("mtkFastAnalyserResult",
		
		contains=c("mtkAnalyserResult")
)

#' The constructor.
#'  @param main a data-frame to hold the main results produced by the Analyser.
#'  @param information a named list to provide supplementary information about the analysis process and its results.

#' @return an object of class \code{\linkS4class{mtkFastAnalyserResult}}
#' @examples mtkmtkFastAnalyserResult()
#' @export mtkmtkFastAnalyserResult
#' @title The constructor

mtkFastAnalyserResult <- function(main, information=NULL) {
	res <- new("mtkFastAnalyserResult", main=main, information=information)
	return(res)
}
