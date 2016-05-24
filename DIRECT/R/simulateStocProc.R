simulateStocProc <-
function (dTime, meanT1, sdT1, meanProc, sdProc, betaProc, MODEL="OU")
{
	L = length (dTime) + 1
	result = rep (0, L)
	result[1] = rnorm (1, mean=meanT1, sd=sdT1)
	
	if (MODEL=="OU")
	{
		sd.tmp = sdProc * (1 - exp (-2*betaProc*dTime)) / 2 / betaProc
		for (l in 2:L)
		{
			mean.tmp = result[l-1] * exp (-betaProc*dTime[l-1]) + meanProc * (1 - exp (-betaProc*dTime[l-1]))
			result[l] = rnorm (1, mean=mean.tmp, sd=sd.tmp[l-1])
		}
	}
	if (MODEL=="RW")
	{
		sd.tmp = sdProc * dTime
		for (l in 2:L)
		{
			mean.tmp = result[l-1] + meanProc * dTime[l-1]
			result[l] = rnorm (1, mean=mean.tmp, sd=sd.tmp[l-1])
		}
	}

	return (result)
}

