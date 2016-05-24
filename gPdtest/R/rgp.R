# This function generates Generalized Pareto random numbers
rgp <- function(n,shape,scale=1)
{
	if(shape!=0)
		return((scale/shape)*(runif(n)**(-shape)-1))
	else return(rexp(n,scale))
}
