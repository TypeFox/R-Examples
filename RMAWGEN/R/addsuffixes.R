NULL


#' 
#' 	Adds suffixes for daily maximum and minimum temperature to the names of a column data frame 
#' 
#' @author Emanuele Cordano, Emanuele Eccel 
#'   
#' @param names a character string vector with column names
#' @param suffix suffixes to add to the first and second groups of column names respectively
#' @param sep separation element
#'  
#'  @export 
#'
#' @seealso \code{\link{getVARmodel}}
#' 
#'       
#' @return  the vector of names with suffixes added 
#' @details This function is used for data frames with duplicated field names
#' 
#' @examples 
#' names <- addsuffixes()




addsuffixes <-
function (names=c("T0001","T0099","T0001","T0099"),suffix=c("_Tx","_Tn"),sep="") {

	ns=length(suffix)
	nn=length(names)
	vn=as.integer(nn/ns)
	
	for (i in 1:ns) names[((i-1)*vn+1):(vn*i)] <- paste(names[((i-1)*vn+1):(vn*i)],suffix[i],sep=sep)
	
	
	
	return(names)
	
}

