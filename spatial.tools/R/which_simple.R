#' Location of Maximum Value
#' 
#' Locates the largest value of the input object.
#' 
#' 
#' @param x a numeric object
#' @param na.rm a logical indicating whether missing values should be removed.
#' @param tie_value A character indicating how to deal with ties.  Can be "NA"
#' (returns an NA if a tie is found) or "random" (returns a single randomly
#' chosen member of the ties if a tie is found) or "first" (returns the first
#' class found).
#' @return An integer of length 1 giving the index of the maximum of x or NA if
#' the maximum of x is not unique, x has no non-NAs, or na.rm=F.
#' @author Jonathan A. Greenberg, Alison R. Mynsberge
#' @seealso \code{\link[base]{which.max}}, \code{\link[base]{which}},
#' \code{\link[base]{max}}
#' @keywords calculate
#' @examples \dontrun{
#' 
#' x<-c(2:4,1,1,NA)
#' y<-c(4,1:3,NA,4)
#' ## The index is only calculated for a unique maximum
#' which.max.simple(x)
#' which.max.simple(y)
#' which.max.simple(y,na.rm=FALSE)
#' which.max.simple(x,na.rm=FALSE)
#' }
which.max.simple=function(x,na.rm=TRUE,tie_value="NA")
{
	if(na.rm)
	{
		x=x[!is.na(x)]
	}
	if(length(x)==0)
	{
		return(NA)
	}
	maxval=max(x)
	if(is.na(maxval))
	{
		return(NA)
	}
	if(sum(x %in% maxval) > 1)
	{
		# Ties exist, figure out what to do with them.
		if(tie_value=="NA")
		{
			return(NA)
		}
		
		if(tie_value=="random")
		{
			tie_postions=which(x==maxval)
			return(sample(tie_postions,size=1))
		}
		
		if(tie_value=="first")
		{
			tie_postions=which(x==maxval)
			return(tie_postions[1])
		}
		
	} else
	{
		return(which.max(x))
	}
}



#' Location of Minimum Value
#' 
#' Locates the smallest value of the input object.
#' 
#' 
#' @param x a numeric object
#' @param na.rm a logical indicating whether missing values should be removed.
#' @param tie_value A character indicating how to deal with ties.  Can be "NA"
#' (returns an NA if a tie is found) or "random" (returns a single randomly
#' chosen member of the ties if a tie is found) or "first" (returns the first
#' class found).
#' @return An integer of length 1 giving the index of the minimum of x or NA if
#' the minimum of x is not unique, x has no non-NAs, or na.rm=F.
#' @author Jonathan A. Greenberg, Alison R. Mynsberge
#' @seealso \code{\link[base]{which.min}}, \code{\link[base]{which}},
#' \code{\link[base]{min}}
#' @keywords calculate
#' @examples \dontrun{
#' 
#' x<-c(4,1:3,NA,4)
#' y<-c(2:4,1,1,NA)
#' ## The index is only calculated for a unique minimum
#' which.min.simple(x)
#' which.min.simple(y)
#' which.min.simple(y,na.rm=FALSE)
#' which.min.simple(x,na.rm=FALSE)
#' }
which.min.simple=function(x,na.rm=TRUE,tie_value="NA")
{
	if(na.rm)
	{
		x=x[!is.na(x)]
	}
	if(length(x)==0)
	{
		return(NA)
	}
	minval=min(x)
	if(is.na(minval))
	{
		return(NA)
	}
	if(sum(x %in% minval) > 1)
	{
		# Ties exist, figure out what to do with them.
		if(tie_value=="NA")
		{
			return(NA)
		}
		
		if(tie_value=="random")
		{
			tie_postions=which(x==minval)
			return(sample(tie_postions,size=1))
		}
		
		if(tie_value=="first")
		{
			tie_postions=which(x==minval)
			return(tie_postions[1])
		}
		
	} else
	{
		return(which.min(x))
	}
}

