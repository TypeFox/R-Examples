logreturns <-
function(x)
{
	log(x[2:(length(x))]/x[1:(length(x)-1)])
}
