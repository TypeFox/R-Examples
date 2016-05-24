# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 8 feb 2011
# licence	: GPL

# Author(s) : Juhui Wang, Hervé Monod, MIA-Jouy en Josas, INRA, 78352
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 300                 $: revision number of the last spread
#' $Author:: jwang            $: author of the last spread
#' $Date:: 2012-08-23 10:26:5#$: date of the last spread

#-----------------------------------------------------------------------

#=======================================================================
# Class: mtkIshigamiEvaluator
#---------------------------------------------------------------------------




#' The Evaluator of Ishigami model defined in the "sensitivity" package
#' @exportClass mtkIshigamiEvaluator
#' @title The mtkIshigamiEvaluator class

setClass(Class="mtkIshigamiEvaluator",
		contains=c("mtkEvaluator"),
		validity = function(object){ object@name=="evaluate"}
)

#' The constructor
#' @param result a data holder to hold the results produced by the process.
#' @return an object of class \code{\linkS4class{mtkIshigamiEvaluator}}
#' @export mtkIshigamiEvaluator

mtkIshigamiEvaluator <- function() {
	res <- new("mtkIshigamiEvaluator",service="Ishigami", ready=TRUE, state=FALSE)
	return(res) 
}



#' Runs the Ishigami Evaluator 
#' @param this an object of class \code{\linkS4class{mtkIshigamiEvaluator}}
#' @param context an object of class \code{\linkS4class{mtkExpWorkflow}}
#' @return invisible()
#' @exportMethod run
#' @title The run method

setMethod(f="run", signature=c(this="mtkIshigamiEvaluator",
				context="mtkExpWorkflow"),
		definition=function(this, context){
			if(this@state) return(invisible())
			
			nameThis <- deparse(substitute(this))
			parameters<-getParameters(this)
			
			X <- context@processesVector$design@result@main
			
			##!! 
			##!! Pre-processing the input data, the processing of the method to implement follows:
			##!!
			out <- ishigami.simule(X,tout=FALSE)
			
			##!!
			##!! post-processing the output of the method:
			##!!
			
			this@result<- mtkEvaluatorResult(main=as.data.frame(out))
			this@state<-TRUE
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		}
)

###################### ORIGINAL CODE ##############################

#' Input factors of the Ishigami model, described in \cite{Saltelli et al., 2000}
#' @docType data
#' @format data.frame of 3 rows (factors) and 4 columns (specifications)
#' @title input factors of the Ishigami model
#' @export ishigami.factors

ishigami.factors <-
		structure(list(nominal = c(0, 0, 0), binf = c(-3.14159265358979, 
								-3.14159265358979, -3.14159265358979), bsup = c(3.14159265358979, 
								3.14159265358979, 3.14159265358979), name = structure(1:3, .Label = c("x1", 
										"x2", "x3"), class = "factor")), .Names = c("nominal", "binf", 
						"bsup", "name"), row.names = c("x1", "x2", "x3"), class = "data.frame")


#' Ishigami model
#' @title Ishigami model, described in Saltelli et al., 2000
#' @param param vector of length 3 or a N x 3 matrix of parameters, with
#'          each parameter ranging between -pi and +pi
#' @return scalar or vector of length N
#' @note calls the function ishigami.fun from the library sensitivity
#' @examples 
#'  ishigami.model( c(-1,0,-1) )
#'  ishigami.model( rbind( c(1,1,1),c(-1,0,-1) )  )
#' @export ishigami.model

ishigami.model <-
		function(param=ishigami.factors$nominal)
{
	if(length(param) == length(ishigami.factors$nominal) ){
		x = matrix(param,nrow=1)
	}
	else x <- param
	return(ishigami.fun(x))
}

#' Simulation of the Ishigami model
#' @title Simulation of the Ishigami model, described in Saltelli et al., 2000
#' @param X matrix or dataframe N x 3 of input values, ranging between -pi and +pi
#' @param tout TRUE to get the input values in the data.frame returned by the function 
#' @return matrice or dataframe if 'tout==TRUE', a vector otherwise
#' @note calls the function ishigami.fun from the library sensitivity
#' @examples 
#'  ishigami.simule( c(-1,0,-1) )
#'  ishigami.simule( rbind( c(1,1,1),c(-1,0,-1) )  )
#' @export ishigami.simule

ishigami.simule <-
		function(X, tout=FALSE)
{
	
	if(is.null(dim(X))) X = matrix(X,nrow=1)
	
	res = t(apply(X,1,ishigami.model))
	res = c(res)  #res=as.data.frame(res)
	if (tout) res  = cbind(X,res)
	return(res)
}
