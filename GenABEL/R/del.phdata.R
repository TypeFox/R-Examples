#' delete phenotypes from phdata
#' 
#' This function is used to delete certain 
#' phenotypes from phenotypic part (phdata) 
#' of an object of \code{\link{gwaa.data-class}}
#' 
#' @param data an object of \code{\link{gwaa.data-class}}
#' @param what which phenotypes (variables) to delete, expressed 
#' as (vector of) names (character) or integer (column of 
#' phdata data frame)
#' @param all if 'all'=TRUE and 'what' is misisng, 
#' all phenotypes are deleted, and only the 'id' 
#' and 'sex' are kept
#' 
#' @author Yurii Aulchenko
#' 
#' @examples
#' data(srdta)
#' phdata(srdta)[1:5,]
#' srdta <- del.phdata(srdta,"qt1")
#' phdata(srdta)[1:5,]
#' srdta <- del.phdata(srdta,all=TRUE)
#' phdata(srdta)[1:5,]
#' 
#' 

del.phdata <- function(data,what,all=FALSE)
{
	if (class(data) != "gwaa.data") stop("data should be of gwaa.data-class");
	if (missing(what) && all) {
		cat("Deleting all data but 'id' and 'sex'\n")
		data@phdata <- data@phdata[,c("id","sex")]
		return(data)
	}
	if (missing(what)) stop("'what' could be missing only if 'all'=TRUE (delete all data)")
	if (is.integer(what)) what <- names(data@phdata)[what]
	
	if (any(what == "id")) stop("you can not delete 'id'")
	if (any(what == "sex")) stop("you can not delete 'sex'")
	
	nms <- names(data@phdata);
	tokeep <- nms[!(nms %in% what)];
	data@phdata <- data@phdata[,tokeep]
	return(data)
}