Linkfunc <-
function(x,rho=0, gamma=0)
{
if(x >= 0 && x < 1) 
temp = 1 - log(1-x) * x^gamma * (1-x)^(rho+1)
else if(x == 1)
temp = 1
else
{
	if(x < 0)
		stop("The first argument in linkfunc() is not negative, so quit")
	else
		stop("The first argument in linkfunc() is larger than 1, so quit!")
}
temp
}

