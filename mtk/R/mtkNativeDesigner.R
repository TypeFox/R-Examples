# Mexico Toolkit
#
# version 	: 0.01
# date		: 27 April 2011
# MAJ   	: 
# licence	: GPL

# Author(s) : Juhui Wang, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 209                 $: revision number of the last spread
#' $Author:: hrichard         $: author of the last spread
#' $Date:: 2011-10-19 17:47:3#$: date of the last spread

#---------------------------------------------------------------------------

##########################################################################
## mtkNativeDesigner class,  definition and methods
###########################################################################

DEBUG.mtkNativeDesigner=FALSE

#' The mtkNativeDesigner class is a sub-class of mtkProcess  used to manage 
#' the sampling task in a simplified way.
#' @title The mtkNativeDesigner class
#' @exportClass mtkNativeDesigner

setClass(Class="mtkNativeDesigner",
		representation=representation(
				design="ANY"
		),
		contains=c("mtkDesigner"),
		validity = function(object){ object@name=="design"}

)

###########################################################################


#' The constructor
#' @param design NULL, a function or a string to specify the designer to use.
#' @param X a data.frame to provide with the results produced off-line.
#' @param information a list to provide with supplementary information about the results produced off-line or
#'  the list of parameters used by the Designer.
#' @return an object of class \code{\linkS4class{mtkNativeDesigner}}
#' @export mtkNativeDesigner

mtkNativeDesigner= function(design=NULL, X=NULL, information=NULL) {
	
	out<-NULL
	st<-FALSE
	if(!is.null(X)) {
		out<- mtkDesignerResult(main=as.data.frame(X), information=information)
		st=TRUE
	}
	res <- new("mtkNativeDesigner",service="Native",design = design,ready=TRUE, state=st, result=out)
	
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
	if(class(design)=="character")
		res<-mtkDesigner(service=design, parameters=p)
	
	return(res) 
}


###########################################################################


#' Generates the design by sampling the factors. 
#' @param this the underlying object of class \code{\linkS4class{mtkNativeDesigner}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
#' @title The run method


setMethod(f="run", signature=c(this="mtkNativeDesigner",context="mtkExpWorkflow"), definition=function(this, context){
			
			if(this@state) return(invisible())
			if(class(this@design)!="function") stop("Please define a sampling function or provide with the off-line sampling results!")
			
			nameThis <- deparse(substitute(this))
			
			res<-this@design 
			this@result <- mtkDesignerResult(main=as.data.frame(res$X), information=res$information)
			this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
			
			
		})
