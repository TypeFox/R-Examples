#' @title Defining vector minimum position  
#' @details
#' Internal function. \code{f.min} is called by \code{splitopt.pls}.
#' @param x vector of value.
#' @param \dots Further arguments passed on to \code{\link{f.min}}. 
#' @return list containing minimun value, the position of the minimun and the values different from NA of a vactor.
#' @keywords internal
#' @export

f.min	<-	function(x,...)
{
	if(!is.null(x))
	{
		v.min	=	min(x[!is.na(x)])
		all.v	=	which(!is.na(x))
		p.min	=	match(v.min,x)
	}
	else
	{
		v.min	=	NULL
		all.v	=	NULL
		p.min	=	NULL
	}	
	list(v.min=v.min,all.v=all.v,p.min=p.min)
}

