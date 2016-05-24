mase = function(forecast, outsampletrue, insampletrue)
{
	if(length(forecast) != length(outsampletrue))
		stop("MASE: the lengths of input vectors must be the same.")
	n = length(insampletrue)
	insamplerr = vector(, (n - 1))
	for(i in 1:(n - 1))
	{
		insamplerr[i] = abs(insampletrue[i+1] - insampletrue[i])
	}
	qt = (outsampletrue - forecast)/(sum(insamplerr)/(n - 1))
	scalederror = mean(abs(qt))
	return(round(scalederror, 6))
}
