calc.rate.lnorm <-
function(par)
{	rate <- 1/exp(par[1] + (par[2]^2)/2)
	return(rate)
}

