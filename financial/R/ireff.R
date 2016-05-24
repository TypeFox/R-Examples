"ireff" <-
function (nom,p) 
{
if (p != Inf) {
	return(((1+nom/(100*p))^p-1)*100)
	}
	else
	{
	return ((exp(nom/100)-1)*100)
	}

}

