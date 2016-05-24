#' Shift bed files
#' 
#' This is a wrapper around the \code{bedcoll} commandline tool.
#' @import rJava
#' @import stats 
#' @importFrom utils write.table read.table
#' @import methods
#' 
#' @param bfile bed filename, without the \code{.bed} extension.
#' @param nshift_min Minimal shift number
#' @param nshift_max Maximal shift number
bedcollr = function(bfile=NULL, nshift_min=1, nshift_max=NULL) {
	bfile = filePath(sprintf("%s.bed", bfile))@path
#	message(sprintf("Shifting bed file: %s", bfile))
	
	paramList = mget(names(formals()),sys.frame(sys.nframe()))
	paramVector = unlist(paramList)
	paramVector = paramVector[!is.null(paramVector)]
	paramVector = str_trim(paramVector)
	
	paramName = names(paramVector)
	names(paramVector) = NULL
	paramName = paste("--", paramName, sep="")
	
	nParam = length(paramName)
	idxOdd = seq(1, nParam * 2, 2)
	idxEven = seq(2, nParam * 2, 2)
	paramNameWithValue = character(nParam * 2)
	paramNameWithValue[idxOdd] = paramName
	paramNameWithValue[idxEven] = paramVector
	system2("bedcoll", paramNameWithValue)
	
	invisible(NULL)
}
