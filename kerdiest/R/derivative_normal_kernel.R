derivative_normal_kernel <-
function(ord, u)
# INPUTS: 
#   "order" order of derivative  
#   "u" array or single value where the kernel is evaluated
{			
    if (ord == 2) # second derivative
		result <- (1/(sqrt(2 * pi)))* exp(-(u^2)/2)*((u^2)-1)
		else if (ord == 4) # fourth derivative
				result <- (1/(sqrt(2 * pi)))* exp(-(u^2)/2) * (3 - (6*(u^2)) + u^4)
		else if (ord == 6) # sixth derivative
				result <- (1/(sqrt(2 * pi)))* exp(-(u^2)/2) * (u^6 - (15*(u^4)) + (45*(u^2)) - 15)
		else if (ord == 8) # eighth derivative
				result <- (1/(sqrt(2 * pi)))* exp(-(u^2)/2) * (u^8 - (28*(u^6)) + (210*(u^4)) - (420*(u^2)) + 105)
		return(result)
	}
