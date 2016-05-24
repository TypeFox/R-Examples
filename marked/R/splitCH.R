#' Split/collapse capture histories
#'
#' splitCH will split a character string vector of capture histories into a matrix. The ch can either be single character or
#' comma separated string. The matrix is appended to the original data set (data) if one is 
#' specified. Will handle character and numeric values in ch. Results will differ depending on content of ch. collapseCH will collapse a
#' capture history matrix back into a character vector. Argument can either be a capture history matrix (chmat) or a dataframe (data)
#' that contains fields with a specified prefix. 
#' 
#' @usage 	splitCH(x="ch", data=NULL, prefix="Time")
#'  
#'	        collapseCH(chmat=NULL, data=NULL, prefix="Time", collapse="")
#' 
#' @aliases splitCH collapseCH
#' @param x A vector containing the character strings of capture histories or the column number or name in the data set \code{data}
#' @param data A data frame containing columnwith value in x if x indicates a column in a data frame
#' @param prefix first portion of field names for split ch
#' @param chmat capture history matrix
#' @param collapse in collapseCH the separator for ch string; defaults to "" but "," also useful if multi-characters are used
#' @return A data frame if data specified and a matrix if vector ch is specified
#' @export splitCH collapseCH
#' @author Devin Johnson; Jeff Laake
#' @examples
#' data(dipper)
#' # following returns a matrix
#' chmat=splitCH(dipper$ch)
#' # following returns the original dataframe with the ch split into columns
#' newdipper=splitCH(data=dipper)
#' # following collapses chmat
#' ch=collapseCH(chmat)
#' # following finds fields in newdipper and creates ch
#' newdipper$ch=NULL
#' newdipper=collapseCH(data=newdipper)

splitCH <- function(x="ch", data=NULL, prefix="Time"){
#   Set value of ch depending on what arguments are set
	if(is.null(data)){
	  ch=x
	} else
	{
		if(!x%in%names(data)) stop(paste("value for data does not contain field", x))
		ch=data[,x]
	}
#   split ch assuming non-numeric fields
	sep=""
	if(length(grep(",",ch[1]))!=0) sep=","
	chmat=do.call("rbind",strsplit(ch,sep))
#   if all fields are numeric split as numeric
	if(!any(!chmat%in%as.character(0:9)))
	chmat=t(sapply(strsplit(ch,sep),function(x)as.numeric(x)))
	if((is.character(x) & length(x)==1) & !is.null(data)){
		colnames(chmat) <- paste(prefix, c(1:ncol(chmat)),sep="")
		rownames(chmat) <- NULL
		return(cbind(data,chmat))
	}
	else{
		colnames(chmat) <- paste(prefix, c(1:ncol(chmat)),sep="")
		return(chmat)
	}
}
collapseCH <- function(chmat=NULL, data=NULL, prefix="Time",collapse=""){
#   Set value of chmat depending on what arguments are set
	if(is.null(data)){
		if(is.null(chmat))
			stop("\nNeither chmat or data were specified.\n")
		else
		if(!is.matrix(chmat))stop("\nchmat must be a matrix\n")
	} else
	{
		fields=names(data)[grep(prefix,names(data))]
		if(length(fields)!=0)
			chmat=subset(data,select=fields)
		else
			stop(paste("\nNo fields found with name containing",prefix,"\n"))
	}
#   collapse ch
	ch=apply(chmat,1,paste,collapse=collapse)
#   if data specified add to dataframe and return; else return character vector
	if(is.null(data))
		return(ch)
	else
	{
		data=cbind(ch=ch,data,stringsAsFactors=FALSE)
		return(data)
	}
}
