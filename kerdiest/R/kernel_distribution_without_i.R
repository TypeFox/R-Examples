kernel_distribution_without_i <-
function(type_kernel,y,x,bw)
# INPUTS:
#   "type" kernel function: "e" Epanechnikov,	"n" Normal, "b" Biweight
#   "y" vector where the kernel estimation is computed
#   "x" sample of data 
#   "bw" bandwidth
# OUTPUT:
# Returns a matrix which stores by columns the estimations for each entry of "y" without one point using the bandwidth stored in "bw" 
{
    n <- length(x)
    AUX <- matrix(0, n, n)
    result <- matrix(0,n,length(y))
	  for(j in 1:length(y))
    { 
      AUX <- matrix(rep.int(outer(y[j],x,"-"),n),nrow=n,byrow=TRUE)
	    aux <- kernel_function_distribution(type_kernel, AUX/bw)
	    diag(aux) <- 0
	    result[,j] <- (1/(n-1))*apply(aux,1,sum)
	  }
  	return(result)
}
