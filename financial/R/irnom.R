"irnom" <-
function (eff,p) 
{
if (p != Inf) {
	return(10^(2*(p-1)/p)*p*(eff+100)^(1/p)-100*p)
	}
	else
	{
	return (100*log((eff+100)/100))
	}

}

