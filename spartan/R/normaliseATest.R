normaliseATest <-
function(result)
{
	if(result<0.5)
	{
		result <- 1-result
	}
	normA <- result
}

