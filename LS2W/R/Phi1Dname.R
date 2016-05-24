`Phi1Dname` <-
function(J, filter.number, family)
{
#
# Program to return a specific character string format
# for a given 1-D discrete father autocorrelation wavelet
#
if(J >= 0.) 
   stop("J must be a negative integer")
return(paste("D1Phi.",  - J, ".", filter.number, ".", family, sep = ""))
}

