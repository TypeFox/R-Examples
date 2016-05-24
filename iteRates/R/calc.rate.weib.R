calc.rate.weib <-
function(par)
{	rate <- par[1]/exp(lgamma(1+1/par[2]))
	return(rate)
}

