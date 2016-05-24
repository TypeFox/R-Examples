gibbs_admkr_nw = function(xh, inicost, k, mutsizp, prob, data_x, data_y)
{
	dm = ncol(data_x)
	data_num = nrow(data_x)
	dv = rn = vector(,dm)
	alpharm = -qnorm(prob/2)
	fx = inicost
	xsum = 0
	for(i in 1:dm)
	{
		rn[i] = rnorm(1)
		xsum = xsum + rn[i]^2
	}
	for(i in 1:dm)
	{
		dv[i] = rn[i]/sqrt(xsum) * rnorm(1) * mutsizp
		xh[i] = xh[i] + dv[i]
	}
	fy = cost_admkr(xh, data_x, data_y)
	r = fx - fy
	if(class(r) != "numeric")
	{
		stop("Please re-run the code.")
	}
	if(r > 0)
	{
        accept = 1
	}
	else
	{
		if(runif(1) < exp(r))
		{
			accept = 1
		}
		else
		{
			accept = 0
		}
	}
	c = mutsizp * ((1 - 1/dm) * sqrt(2 * pi) * exp(alpharm^2/2))/(2 * alpharm) + 1/(dm * prob * (1 - prob))
	if(accept == 1)
	{
		accept_h = 1
		inicost  = fy
		mutsizp = mutsizp + c * (1 - prob)/k
	}
	else
	{
		accept_h = 0
		for(i in 1:dm)
		{
			xh[i] = xh[i] - dv[i]
		}
		mutsizp = mutsizp - c * prob/k
	}
	return(list(x = xh, cost = inicost, accept_h = accept_h, mutsizp = mutsizp))
}
