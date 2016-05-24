rmsse = function(forecast, outsampletrue, insampletrue)
{
	if(length(forecast) != length(outsampletrue))
		stop("RMSSE: the lengths of input vectors must be the same.")
	n = length(insampletrue)
	insamplerr = vector(, (n - 1))
	for(i in 1:(n - 1))
	{
		insamplerr[i] = abs(insampletrue[i+1] - insampletrue[i])
	}
	qt = (outsampletrue - forecast)/(sum(insamplerr)/(n - 1))
	scalederror = sqrt(mean(qt^2))
	return(round(scalederror, 6))
}

