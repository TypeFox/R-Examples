`Psi1Dname` <-
function(J, filter.number, family)
{
#
# Program to return a specific character string format
# for a given 1-D Discrete autocorrelation wavelet
#
if(J >= 0.) 
   stop("J must be a negative integer")
return(paste("D1Psi.",  - J, ".", filter.number, ".", family, sep = ""))
}

