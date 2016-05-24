`estnoise` <-
function(y, yh, regression = FALSE, nmse = TRUE)
{
	# test if y and yh have same size
	if(length(y) != length(yh))
	{
		print("Lengths of label vector y and denoised label vector yh must be equal")
	}
	# test if y and yh are column vectors and convert them if necessary
	if(nrow(y) != length(y))
	{
		y <- matrix(y, length(y), 1)
	}
	if(nrow(yh) != length(yh))
	{
		yh <- matrix(yh, length(yh), 1)
	}
	
	n <- length(y)
	
	if(regression == TRUE)
	{
		# mean squared error
		mse <- sum((y - yh)^2)
		if(nmse == FALSE)
		{
			# conventional mean squared error
			return((1/n)*mse)
		}
		else
		{
			# normalized mean squared error
			return(mse/sum( (y - (1/n)*sum(y))^2 ))
		}
	}
	else
	{
		# 0-1 loss
		return((1/n)*sum(y != yh))
	}
}

