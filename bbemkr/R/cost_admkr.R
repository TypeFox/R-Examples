cost_admkr = function(x, data_x, data_y)
{
	data_num = nrow(data_x)
	dm = ncol(data_x)
	hn = data_num^(-1/(4 + dm))
	bn = data_num^(-3/(2 * dm + 11))
	b = exp(0.5 * x[dm+1]) * bn
	data_e = residu(x, data_x = data_x, data_y = data_y)
    cont = exp(-0.5 * log(2.0 * pi))
    tau2 = vector(, dm + 1)
    for(i in 1:(dm+1))
    {
        tau2[i] = exp(x[i])
    }
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
	for(i in 2:(data_num - 1))
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
	# log Jacobi
	for(i in 1:(dm+1))
	{
		logf = logf + x[i]
	}
	# log priors
	logf = logf + logpriorh2(tau2[dm + 1] * bn^2) + log(bn^2)
	for(i in 1:dm)
	{
		logf = logf + logpriorh2(tau2[i] * hn^2) + log(hn^2)
	}
	return(-logf)
}
