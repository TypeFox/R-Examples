ccChecks <-
function(nCC,
										 threshold=NULL)
{
	##
	if(min(nCC) < 0 & length(nCC) == 1)
		return("* Case-control sample size 'nCC' is negative")
	if(min(nCC) < 0 & length(nCC) > 1)
		return("* At least one case-control sample size 'nCC' is negative")
	
	##
	if(!is.null(threshold))
	{
		if(length(threshold) != 2)
			return("* 'threshold' is not a pair of numbers")
	}
  
	##
	return("")
}
