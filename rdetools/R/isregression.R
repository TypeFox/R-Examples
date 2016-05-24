`isregression` <-
function(y, regression = FALSE)
{
	if(regression != TRUE)
	{
		# find out whether this is a regression or classification problem
		if(sum(y == -1) + sum(y == 1) < length(y))
		{
			# there are labels which are not element of {-1, 1} => regression problem
			regression <- TRUE
		}
	}
	
	return(regression)
}

