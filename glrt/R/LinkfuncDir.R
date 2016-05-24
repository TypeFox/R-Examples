LinkfuncDir <-
function(x,rho=0, gamma=0)
{
if(x >= 0 && x < 1) 
temp = (1 - x)^rho * x^gamma + log(1-x) * ( (rho + 1) * (1-x)^rho * x^gamma - gamma * x^(gamma-1) * (1 - x)^(rho+1) )
else if(x == 1)
temp = 1
else
{
warning("first argument in linkfunc() is not in [0, 1]!")
return(NULL)
}
temp
}

