# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 8 feb 2011
# licence	: GPL
# Author(s) : Juhui WANG, based on a version of H. Richard and H. MONOD
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 275                 $: revision number of the last spread
#' $Author:: hrichard         $: author of the last spread
#' $Date:: 2012-01-17 11:31:0#$: date of the last spread

#---------------------------------------------------------------------------
# Nota : 
# fusion r29 et r24 (r24 est le bon code)
# 
# AllowedType is a global variable, cadidate to be initialized by a
# form of preprocessing (SEE IT LATER)
# allowedTypes : character, double, integer, logical
#===================================================================

# class Declaration

#' mtkValue is a class to manage the mapping between the value as a string  of a variable  and  its value utilisable by the programming language.
#' @slot name name of the variable.
#' @slot type data type of the variable.
#' @slot val value of the variable in the right type (i.e. character, double or logical).
#' @exportClass mtkValue vector
#' @title The mtkValue class

setClass(Class= "mtkValue", representation=representation(
                              name="character",
                              type="character",
                              val="ANY"),
         prototype=prototype(name="unknown",type='', val=NULL), 

         validity=function(object){
           goodType <- switch(object@type,
                              character = TRUE,
							  numeric = TRUE,
                              double = TRUE,
                              integer = TRUE,
                              logical = TRUE,
							  string = TRUE,
							  float = TRUE,
							  null = TRUE,
							  NULL = true,
                              FALSE
                              )
							  
           if(!goodType)
             {
               stop("Wrong type! Only following types are allowed: character, string, double, float, integer, logical")
             }
           return(TRUE)
         }
         ##         contains=character(),
         )


#' The constructor method 
#' @return an object of class \code{\linkS4class{mtkParameter}}
#' @export mtkParameter

mtkValue = function(name='unknown', type='', val=NULL) {
	if(type == '') type <- typeof(val)
	if(type=='NULL') type <- 'null'
	if(type == 'string') type <- 'character'
	if(type == 'float') type <- 'double'
	cmd <- paste("as.", type, sep="")
	val <- eval(do.call(cmd, list(val)))
	res <- new("mtkValue", name=name, type=type, val=val)
	return(res)
}


#===================================================================
# Method definition
#' The getName method
#' @param this the underlying object of class \code{\linkS4class{mtkValue}}
#' @return the name of the variable.
#' @exportMethod getName

# return the "name"  of the object mtkValue
# move in mtkAllGenerics.R, uncomment otherwise :
# setGeneric ("getName",   function(this) {standardGeneric("getName")} )

setMethod(f="getName",
		signature=c("mtkValue"),
		definition=function(this) {
			return(this@name)
		}
)

#===================================================================

#' The getValue method
#' @param this the underlying object of class \code{\linkS4class{mtkValue}}
#' @return the value of the variable.
#' @exportMethod getValue
# return the real value of the mtkValue class object

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getValue", function(this) {standardGeneric("getValue")} )

# retourne une valeur associée à son nom

setMethod(f="getValue",
		signature=c("mtkValue"),
		definition=function(this) {
			cmd <- paste("as.", this@type, sep="")
			val <- eval(do.call(cmd, list(this@val)))
			return (val)
		})
#===================================================================
# Method definition
#' The getValueType method
#' @param this the underlying object of class \code{\linkS4class{mtkValue}}
#' @return the data type of the variable.
#' @exportMethod getValueType

### return the value type of the mtkValue class object

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getValueType", function(this) {standardGeneric("getValueType")} )

setMethod(f="getType",
		signature=c("mtkValue"),
		definition=function(this) {
			return(this@type)
		})

#===================================================================
# Method definition
#' The getValueType method
#' @param this the underlying object of class \code{\linkS4class{mtkValue}}
#' @return the data type of the variable.
#' @exportMethod getValueType

### return the value type of the mtkValue class object

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getValueType", function(this) {standardGeneric("getValueType")} )

setMethod(f="setName", signature=c(this="mtkValue", name="character"),
		definition=function(this, name ) {
			nameThis <- deparse(substitute(this))
			this@name <- name
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

#===================================================================
# Method definition
#' The getValueType method
#' @param this the underlying object of class \code{\linkS4class{mtkValue}}
#' @return the data type of the variable.
#' @exportMethod getValueType

### return the value type of the mtkValue class object

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getValueType", function(this) {standardGeneric("getValueType")} )

setMethod(f="setValue", signature=c(this="mtkValue", val="ANY"),
		definition=function(this, val ) {
			nameThis <- deparse(substitute(this))	
			
			cmd <- paste("as.", this@type, sep="")
			this@val <- eval(do.call(cmd, list(val)))
			
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})
#===================================================================
# Method definition
#' The getValueType method
#' @param this the underlying object of class \code{\linkS4class{mtkValue}}
#' @return the data type of the variable.
#' @exportMethod getValueType

### return the value type of the mtkValue class object

# move in mtkAllGenerics.R, uncomment otherwise :
#setGeneric ("getValueType", function(this) {standardGeneric("getValueType")} )

setMethod(f="setType", signature=c(this="mtkValue", type="character"),
				definition=function(this, type ) {
			nameThis <- deparse(substitute(this))	
			if(type == '') type <- typeof(this@val)
			if(type=='NULL') type <- 'null'
			if(type == 'string') type <- 'character'
			if(type == 'float') type <- 'double'
			this@type <-  type
			cmd <- paste("as.", type, sep="")
			this@val <- eval(do.call(cmd, list(this@val)))
			assign(nameThis, this, envir=parent.frame())
			return(invisible())
		})

##########################################################################
## Juhui WANG le 28/08/2012
#' Prints the information managed by the class.
#' @param x the underlying object of class \code{\linkS4class{mtkValue}}.
#' @return a summary of the data managed by the \code{\linkS4class{mtkValue}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="print", signature="mtkValue", definition=function(x,...){
			cat("---------------------------------------------------------------------------\n")
			cat(paste("	 HERE IS INFORMATION ABOUT THE OBJECT OF mtkValue ",": \n", sep=""))
			cat("---------------------------------------------------------------------------\n\n\n")
			cat("Name : ", x@name, "\n")
			cat("Type : ", x@type, "\n")
			cat("Value : ", unlist(x@val), "\n\n")
			invisible()
		}
)

##########################################################################
## Juhui WANG le 28/08/2012
#' Prints the information managed by the class.
#' @param x the underlying object of class \code{\linkS4class{mtkValue}}.
#' @return a summary of the data managed by the \code{\linkS4class{mtkValue}}.
#' @exportMethod  print
#' @title The print method
setMethod(f="show", signature="mtkValue", definition=function(object){
			print(object)
			invisible()
		}
)

