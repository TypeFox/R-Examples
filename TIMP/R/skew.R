"skew" <- 
function(numax, deltanu, b, nu, nupower=1)
{
	arg <- 1 + (2 * b * (nu - numax))/deltanu
        res <- arg
        res[which(arg>0)] <- exp( - log(2) * (log(arg[which(arg>0)])/b)^2)
        res[which(arg<=0)] <- 0 
	if(nupower!=1)
		res <- res * nu^nupower 

	res
}
