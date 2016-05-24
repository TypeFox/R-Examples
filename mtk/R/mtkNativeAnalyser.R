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
#' $Rev:: 223                 $: revision number of the last spread
#' $Author:: hmonod           $: author of the last spread
#' $Date:: 2011-11-04 00:25:0#$: date of the last spread

#---------------------------------------------------------------------------

##########################################################################
## mtkNativeAnalyser class,  definition and methods
###########################################################################

DEBUG.mtkNativeAnalyser=FALSE


#' The mtkNativeAnalyser class is a sub-class of mtkProcess  used to manage 
#' the analysis task for the cases where either the results are produced off-line or
#' the Analyser is user-defined with the help of a R function.
#' @title The mtkNativeAnalyser class
#' @exportClass mtkNativeAnalyser

setClass(Class="mtkNativeAnalyser",
		representation=representation(
				analyze="ANY"
		),
		contains=c("mtkAnalyser"),
		validity = function(object){ object@name=="analyze"}

)

###########################################################################


#' The constructor
#' @param Analyser NULL, a function or a string to specify the Analyser to use.
#' @param X a data-frame to provide with the results produced off-line.
#' @param information a named list to provide with supplementary information about the results produced off-line or the parameters 
#' used by the Analyser.
#' @return an object of class \code{\linkS4class{mtkNativeAnalyser}}
#' @export mtkNativeAnalyser

mtkNativeAnalyser= function(analyze=NULL, X=NULL, information=NULL) {
	
	out<-NULL
	st<-FALSE
	if(!is.null(X)) {
		out<- mtkAnalyserResult(main=as.data.frame(X), information=information)
		st=TRUE
	}
	res <- new("mtkNativeAnalyser",service="Native",analyze=analyze,ready=TRUE, state=st, result=out)
	p<-NULL
	if(!is.null(information)){
### HM, 5/6/11: debut
###		p<-sapply(1:length(information), function(i,x,lab)
###				{ 
###					mtkParameter(valName=lab[i],type=class(x[[i]]), valString=as.character(x[[i]]))
###				}, information, names(information)
###		)
          p <- make.mtkParameterList(information)
        }
### HM, 5/6/11: fin
	if(class(analyze)=="character")
		res <- mtkAnalyser(service=analyze, parameters=p)
	
	return(res) 
}


###########################################################################


#' Runs the sensitivity analysis. 
#' @param this the underlying object of class \code{\linkS4class{mtkNativeAnalyser}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
#' @title The run method


setMethod(f="run", signature=c(this="mtkNativeAnalyser",context="mtkExpWorkflow"), definition=function(this, context){
			
			if(this@state) return(invisible())
			
			if(class(this@analyze)!="function") stop("Without the results produced off-line, you must provide with a analysis method!")
								
			nameThis <- deparse(substitute(this))
			
			res<-this@analyze 
			this@result <- mtkAnalyserResult(main=as.data.frame(res$X), information=res$information)
			this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
			
			
		})
