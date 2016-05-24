# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 8 feb 2011
# licence	: GPL

# Author(s) : Juhui WANG, 10 April, 2014
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/

#=======================================================================
#' The mtkLevels class make.mtkFeature
#' @slot levelValue value of the parameter among character, double, integer, logical
#' @slot weight numeric optional wheighting (bounded between 0 to 1)
#' @exportClass mtkLevels vector
#===================================================================
# class Declaration

setClass(Class="mtkLevels",
		representation=representation(
				type = "character",
				levels="vector",
				weights="numeric"
		),

	prototype=prototype(type="categorical", levels=vector(), weights =numeric(0))
		
)

#===================================================================
#' The constructor method
#' @param levelValue a character, double, integer, or logical atomic object
#' @param weight a numeric value, default weight=1 means no weight
#' @return an object of class \code{\linkS4class{mtkLevel}}
#' @export mtkLevel
mtkLevels = function(type = "categorical", levels=vector(), weights=numeric(0))
{
	l <- length(levels)
	w <- length(weights)
	if(w==0 && l> 0) weights <- rep(1/l,l)
	if(w != 0 && w != l) stop("Every levels must be weighted!\n")
	
	new(Class="mtkLevels", type=type, levels=levels, weights=weights)
}

###########################################################################
#' Gives a name to the type
#' @param this the underlying object of class \code{\linkS4class{mtkLevels}}
#' @param name a string indicating the processing step associated with this process. 
#' @return invisble()
#' @exportMethod setName
#' @title The setType method

setMethod(f="setType", signature=c(this="mtkLevels", type="character"),
		definition=function(this, type ) {
			nameThis <- deparse(substitute(this))
			this@type <-  type
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

###########################################################################
#' get the type
#' @param this the underlying object of class \code{\linkS4class{mtkLevels}}
#' @return a string indicating the processing step associated with this mtkLevels.
#' @exportMethod getType
#' @title The getType method

setMethod(f="getType", signature=c(this="mtkLevels"),
		definition=function(this) {
			return(this@type)
		})


###########################################################################
#' Gives a new vector of levels
#' @param this the underlying object of class \code{\linkS4class{mtkLevels}}
#' @param name a string indicating the processing step associated with this process. 
#' @return invisble()
#' @exportMethod setName
#' @title The setType method

setMethod(f="setLevels", signature=c(this="mtkLevels", levels="vector"),
		definition=function(this, levels) {
			nameThis <- deparse(substitute(this))
			l <- length(levels)
			w <- length(this@weights)
			if(l != w) cat(" Warning: the numbers of  'levels' and 'weights' are not equal! \n Don't forget to modify the 'weights' too.\n" )
			this@levels <-  levels
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

###########################################################################
#' get the levels
#' @param this the underlying object of class \code{\linkS4class{mtkLevels}}
#' @return a string indicating the processing step associated with this mtkLevels.
#' @exportMethod getLevels
#' @title The getLevels method

setMethod(f="getLevels", signature=c(this="mtkLevels"),
		definition=function(this) {
			return(this@levels)
		})

###########################################################################
#' Gives a new vector of weights
#' @param this the underlying object of class \code{\linkS4class{mtkLevels}}
#' @param name a string indicating the processing step associated with this process. 
#' @return invisble()
#' @exportMethod setName
#' @title The setType method

setMethod(f="setWeights", signature=c(this="mtkLevels", weights="numeric"),
		definition=function(this, weights ) {
			nameThis <- deparse(substitute(this))
			l <- length(this@levels)
			w <- length(weights)
			if(l != w) stop(" The numbers of  'levels' and 'weights' are not equal! \n Modify first the 'levels'  and then the 'weights'.\n" )
			this@weights <-  weights
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

###########################################################################
#' get the weights
#' @param this the underlying object of class \code{\linkS4class{mtkLevels}}
#' @return a string indicating the processing step associated with this mtkLevels.
#' @exportMethod getWeights
#' @title The getweights method

setMethod(f="getWeights", signature=c(this="mtkLevels"),
		definition=function(this) {
			return(this@weights)
		})

##########################################################################
## Juhui WANG le 28/08/2012
#' Prints the information managed by the class.
#' @param x the underlying object of class \code{\linkS4class{mtkLevels}}.
#' @return a summary of the data managed by the \code{\linkS4class{mtkLevels}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="print", signature=c(x="mtkLevels"), definition=function(x,...){
			
			levels <- unlist(x@levels)
			weights <- unlist(x@weights)
			
			tmp <- data.frame(levels, weights)
			if(length(levels)>0){
				cat("type: ",x@type,"\n")
				print(tmp,...)
			}
			return(invisible())
		}
)

##########################################################################
## Juhui WANG le 28/08/2012
#' Prints the information managed by the class.
#' @param x the underlying object of class \code{\linkS4class{mtkLevels}}.
#' @return a summary of the data managed by the \code{\linkS4class{mtkLevels}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="summary", signature=c(object="mtkLevels"), definition=function(object,...){
			
			cat("---------------------------------------------------------------------------\n")
			
			cat(paste("	 HERE IS A SUMMARY OF THE  \"", object@name, "\": \n", sep=""))
			show(object)
		}
)

# Method definition
#' The mtkLevel show method
#' @param this the underlying object of class \code{\linkS4class{mtkLevel}}
#' @return invisible.
#' @exportMethod show

### return nothing but the write on the R console the content of the mtkLevel

setMethod(f="show",
		signature=c("mtkLevels"),
		definition=function(object) {
			
				print(object)
			
		})


## )
#===================================================================
