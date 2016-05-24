loglikelihood_admkr = function(h2, data_x, data_y)
{
    dm = ncol(data_x)
    data_num = nrow(data_x)
    bn = data_num^(-3/(2 * dm + 11))
    b = sqrt(h2[dm+1]) * bn
    xp = vector(,(dm+1))
    for(i in 1:(dm+1))
    {
        xp[i] = log(h2[i])
    }
    data_e = residu(xp, data_x, data_y)
    cont = exp(-0.5 * log(2.0 * pi))
	epsilon = scale(data_e)
	std = sd(data_e)
	nsum = xsum = 0
	for(j in 2:data_num)
	{
		temp = epsilon[1] - epsilon[j]
		if(abs(temp) > 0)
		{
			xsum = xsum + cont * exp(-0.5 * (temp/b)^2)/b
			nsum = nsum + 1
		}
	}
	logf = log(xsum/nsum/std)
	for(i in 2:(data_num-1))
	{
		nsum = xsum = 0
		for(j in 1:(i-1))
		{
			temp = epsilon[i] - epsilon[j]
			if(abs(temp) > 0)
			{
				xsum = xsum + cont * exp(-0.5 * (temp/b)^2)/b
				nsum = nsum + 1
			}
		}
		for(j in (i+1):data_num)
		{
			temp = epsilon[i] - epsilon[j]
			if(abs(temp) > 0)
			{
				xsum = xsum + cont * exp(-0.5 * (temp/b)^2)/b
				nsum = nsum + 1
			}
		}
		logf = logf + log(xsum/nsum/std)
	}
	nsum = xsum = 0
	for(j in 1:(data_num-1))
	{
		temp = epsilon[data_num] - epsilon[j]
		if(abs(temp) > 0)
		{
			xsum = xsum + cont * exp(-0.5 * (temp/b)^2)/b
			nsum = nsum + 1
		}
	}
	logf = logf + log(xsum/nsum/std)
    return(logf)
}











