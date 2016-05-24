symmetric.alpha.divergence <-
function(x, y)
{
	warning("Call to symmetric.alpha.divergence(...) is deprecated!");
	return(symmetricAlphaDivergence(x,y));
}


symmetricAlphaDivergence <-
function(x, y)
{
	return ( 4*(1-sum(sqrt(x*y)) ))
	
}
