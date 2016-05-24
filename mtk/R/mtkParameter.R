# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 8 feb 2011
# licence	: GPL

# Author(s) : Juhui WANG, based on a version of H. Richard INRA MIA BioSP
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 223                 $: revision number of the last spread
#' $Author:: hmonod           $: author of the last spread
#' $Date:: 2011-11-04 00:25:0#$: date of the last spread
#
#' Nota :
#' fusion r28 et r29 (code de la r28) 
#---------------------------------------------------------------------------
#' inherit of : mtkValue
#' @exportClass mtkParameter
#' @title The mtkParameter class
setClass(Class = "mtkParameter",
		representation = representation(),
		contains = c("mtkValue")
)
#===================================================================
#' The constructor method 
#' @param ... an argument for building a mtkValue (see mtkValue Constructor)
#' @return an object of class \code{\linkS4class{mtkParameter}}
#' @export mtkParameter

mtkParameter = function(name='unknown', type='logical', val=NULL) {
		if(type == '') type <- typeof(val)
		if(type=='NULL') type <- 'null'
		if(type == 'string') type <- 'character'
		if(type == 'float') type <- 'double'
		cmd <- paste("as.", type, sep="")
		val <- eval(do.call(cmd, list(val)))
		res <- new("mtkParameter", name=name, type=type, val=val )
		return(res)
	}




# #===================================================================
## HM, 5/6/11: internal function
## Makes a list of mtkParameter elements from a simple named list
## HR 2011-08-22 convert function in method make.mtkParameterList
## HM 2011-11-03 conservation du renommage, retour a une fonction
#===================================================================
# 
#' Makes a list of mtkParameter elements from a simple named list
#' @param list a named list
#' @return a list of mtkParameter objects
#' @examples make.mtkParameterList(list(min=-1,max=+1,shape="hello"))
## setMethod(f="make.mtkParameterList", 
## 		signature=c("list"),
## 		definition= function(list=list()) {
make.mtkParameterList <- function(x=list()) {
  nbElements <- length(x)
  listMtkElements <-list()
  tmpMtkElement <- NULL
  
  for (i in 1:nbElements) {
    #nameElement <-names(x[i])
    #valElement <- x[[i]]
    #assign(nameElement, valElement, pos=1)
    #tmpMtkElement <- mtkParameter(nameElement)
    tmpMtkElement <- mtkParameter()
    tmpMtkElement@name <- names(x[i])
    tmpMtkElement@val <- x[[i]]
    tmpMtkElement@type <- typeof(x[[i]])
    
    listMtkElements <- c(listMtkElements, tmpMtkElement)
  }
  
  return(listMtkElements)
}

