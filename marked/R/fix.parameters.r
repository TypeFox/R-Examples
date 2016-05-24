#' Fixing real parameters in crm models
#' 
#' Creates matrix with appropriate structure for fixing real parameters in a
#' crm model for a specified set of capture histories and occasions.
#' 
#' 
#' @param x processed data list
#' @param subset logical expression for set of animals; used in subset function
#' @param occasions vector of occasion numbers
#' @param values either single or vector of values; if latter, must match
#' length of occasions
#' @return Matrix with 3 columns: 1) animal row number, 2) occasion, 3) value
#' for real parameter
#' @author Jeff Laake
#' @export
fix.parameters=function(x,subset=NULL,occasions,values)
#   Creates matrix for fixing real parameters; matrix has 3 columns: 
#	1) animal row number, 2) occasion, 3) value for real parameter
#   This function creates that matrix for a specified subset of animals (subset),
#   vector of occasions (occasions) and sets to specific values (values).
#
#   Arguments:
#    x         : processed data list
#    subset    : logical expression for set of animals; used in subset function
#    occasions : vector of occasion numbers
#    values    : either single or vector of values; if latter, must match length of occasions	
# 
#   Value:  matrix with 3 columns as described above and a row for each selected animal-occasion.
{
	x$xid=1:nrow(x$data)
	if(is.null(subset))
		ids=1:nrow(x$data)
	else
	    ids=subset(x$data,subset=subset,select="xid")
	if(length(values)==1)
	   fixed=cbind(rep(ids,length(occasions)),rep(occasions,each=length(ids)),rep(values,length(ids)*length(occasions)))
    else
	   if(length(values)==length(occasions))
	      fixed=cbind(rep(ids,length(occasions)),rep(occasions,each=length(ids)),rep(values,each=length(ids)))
	   else
		  stop("Error in fix.parameters: Number of values specified does not match number of occasions")
	return(fixed)
}
