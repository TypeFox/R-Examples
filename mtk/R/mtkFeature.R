# Mexico Toolkit
#
# version 	: 0.01
# date		: 30 nov 2009
# MAJ   	: 14 janv 2010
# licence	: GPL

# Author(s) : Hervé Richard MIA INRA BioSP
# repository: https://mulcyber.toulouse.inra.fr/projects/baomexico/
# Web		: http://www.reseau-mexico.fr/
#
#' $Rev:: 295                 $: revision number of the last spread
#' $Author:: hmonod           $: author of the last spread
#' $Date:: 2012-07-17 10:29:1#$: date of the last spread

#-----------------------------------------------------------------------
#' Nota :
#' fusion rev 26 et rev 29 (R26 contenant les bonnes modif)
#=======================================================================
#' The mtkFeature class
#' inherit of : mtkValue
#' @exportClass mtkFeature
#' 
setClass(Class = "mtkFeature",
#		representation = representation(),
		contains = c("mtkValue")
)
#===================================================================
#' The constructor
#' @param ... an argument for building a mtkValue (see mtkValue Constructor)
#' @return an object of class \code{\linkS4class{mtkFeature}}
#' @export mtkFeature
#' 
mtkFeature <- function(name='unknown', type='logical', val=NULL)
{	if(type == '') type <- typeof(val)
	if(type == 'NULL') type <- 'null'
	cmd <- paste("as.", type, sep="")
	val <- eval(do.call(cmd, list(val)))
	res <- new("mtkFeature", name=name, type=type, val=val)
	return(res)
	
}
#===================================================================


#===================================================================
### HM, 5/6/11: internal function
## HR 2011-10-06 convert function in method
## HM 2011-11-03 conservation du renommage, retour a une fonction
#===================================================================
# 
#' Makes a list of mtkFeature elements from a simple named list
#' @param list a named list
#' @return a list of mtkFeature objects
#' @examples make.mtkFeatureList(list(min=-1,max=+1,shape="hello"))


make.mtkFeatureList <- function(x=list()) {
	nbElements <- length(x)
	listMtkElements <-list()
	tmpMtkElement <- NULL
	
	for (i in 1:nbElements){
		#nameElement <-names(list[i])
		#valElement <- list[[i]]
		#assign(nameElement, valElement, pos=1)
		#tmpMtkElement <- mtkFeature(nameElement)
		tmpMtkElement <- mtkFeature()
		tmpMtkElement@name <- names(x[i])
		tmpMtkElement@val <- x[[i]]
		tmpMtkElement@type <- typeof(x[[i]])
		
		listMtkElements <- c(listMtkElements, tmpMtkElement)
	}
	return(listMtkElements)
}


