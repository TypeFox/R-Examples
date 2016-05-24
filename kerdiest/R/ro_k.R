ro_k <-
function(type_kernel)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,	"n" Normal, 
#                                  "b" Biweight, "t" Triweight         

{
	if(type_kernel == "e")	
    result <- 2*0.12857	
	else 
		if(type_kernel == "n")	
      result <- 2*0.28209		 
	else 
		if(type_kernel == "b")	
		  result <- 2*0.10823	 
	else 
		if(type_kernel == "t")	
      result <- 2*0.095183	
  return(result)
}
