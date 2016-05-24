kernel_function_distribution <-
function(type_kernel,u)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,	"n" Normal, 
#                                  "b" Biweight, "t" Triweight         
#   "u" array or single value where the kernel is evaluated

{
	if(type_kernel == "e")	
	{
	  result <- u
        Logic0 <- (u <= -1)
        Logic1 <- (u >= 1)
        Logic2 <- (u > -1 & u < 1)
        result[Logic0] <- 0
        result[Logic1] <- 1
	  Uval <- result[Logic2]
    result[Logic2] <- 0.75 * Uval * (1 - (Uval^2)/3) + 0.5
    return(result)
	}
	
	else 
		if(type_kernel == "n")	
		{
		result <- pnorm(u)
		return(result)
		}
		 
	else 
		if(type_kernel == "b")	
		{
		  result <- u
         	Logic0 <- (u <= -1)
          Logic1 <- (u >= 1)
       	  Logic2 <- (u > -1 & u < 1)
          result[Logic0] <- 0
	        result[Logic1] <- 1
		  Uval <- result[Logic2]
      result[Logic2] <- (((15/16) * Uval) -((5/8)* (Uval^3))+((3/16)* (Uval^5))) + 0.5
      return(result)
	  }
		 
	else 
		if(type_kernel == "t")	
		{
		  result <- u
         	Logic0 <- (u <= -1)
          Logic1 <- (u >= 1)
       	  Logic2 <- (u > -1 & u < 1)
       	  Logic3 <- (u > -1 & u < 1)
          result[Logic0] <- 0
	        result[Logic1] <- 1
		  Uval <- result[Logic2]
      result[Logic2] <- (((35/32) * Uval) -((35/32)* (Uval^3))+((21/32)* (Uval^5)))-((5/32)*(Uval^7)) + 0.5
      return(result)
	  }

}
