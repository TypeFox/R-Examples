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
#' $Rev:: 296                 $: revision number of the last spread
#' $Author:: hmonod           $: author of the last spread
#' $Date:: 2012-07-17 10:30:3#$: date of the last spread

#---------------------------------------------------------------------------

DEBUG.mtkNativeEvaluator=FALSE



#' The mtkNativeEvaluator class is a sub-class of the class \code{\linkS4class{mtkEvaluator}} which makes the interactive simulation task easier.
#' @title The mtkNativeEvaluator class
#' @exportClass mtkNativeEvaluator

setClass(Class="mtkNativeEvaluator",
		representation=representation(
				model="ANY"
		),
		contains=c("mtkEvaluator"),
		validity = function(object){ object@name=="evaluate"}

)

###########################################################################
## Methods definition
###########################################################################
#' @param model NULL or a function to specify the model to simulate.
#' @param Y a data-frame to provide the off-line produced simulation.
#' @param information a named list to give supplementary information about the process (parameters, attributes, etc.).
#' @return an object of class \code{\linkS4class{mtkNativeEvaluator}}
#' @export mtkNativeEvaluator

mtkNativeEvaluator= function(model=NULL, Y=NULL, information=NULL) {
	
	re<-NULL
	st<-FALSE
	if(!is.null(Y)) {
		re<- mtkEvaluatorResult(main=Y, information=information)
		st<-TRUE
	}
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
	## if(class(model)=="character"|is.function(model)) ## modified by HM, 5/2/12
	res <- new("mtkNativeEvaluator", service="Native", model=model, parameters=p, state=st, result=re)
	
	
	if(class(model)=="character") ## refait par JW 03/08/2013
		res<-mtkEvaluator(service=model, parameters=p)
	return(res) 
}


###########################################################################

#' Launches the model simulation. 
#' @param this the underlying object of class \code{\linkS4class{mtkNativeEvaluator}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
#' @title The run method


setMethod(f="run", signature=c(this="mtkNativeEvaluator",context="mtkExpWorkflow"), definition=function(this, context){
			
			if(this@state) return(invisible())
			
			if(class(this@model)!="function") stop("Please define a model or provide with the results of off-line simulation")
			
			nameThis <- deparse(substitute(this))
			
			pstr <- NULL
			p <- getParameters(this)
			if(!is.null(p)) pstr <- paste(" ,",paste(names(p), p, sep='=',collapse=','))
			
			X <- extractData(context, name=c('design'))
			out <- NULL
			
			for(i in 1:nrow(X)){
				D <- X[i,]
				func <- paste("this@model","(", paste(names(D),D,sep='=',collapse=","), 
						pstr, ")", sep="")
			
				out <- c(out,eval(parse(text=func)))
			}
				
			res <- list(Y=out,information=list(factors=colnames(X)))
			
			this@result <- mtkEvaluatorResult(main=as.data.frame(res$Y), information=res$information)
			this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
			
			
		})
