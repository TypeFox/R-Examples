A1_k <-
function(type_kernel)
# INPUTS:
#   "type_kernel" kernel function: "e" Epanechnikov,	"n" Normal, 
#                                  "b" Biweight, "t" Triweight         

{
	if(type_kernel == "e")	
    result <- 0.12857	
	else 
		if(type_kernel == "n")	
      result <- 0.28209		 
	else 
		if(type_kernel == "b")	
		  result <- 0.10823	 
	else 
		if(type_kernel == "t")	
      result <- 0.095183	
  return(result)
}
