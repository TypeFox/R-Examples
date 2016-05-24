`Psi2Dname` <-
function(J, filter.number, family, switch)
{
#
# Program to return a specific character string format
# for a given 2-D Discrete autocorrelation wavelet
#
if(J >= 0.)
   stop("J must be a negative integer")
if(switch == "direction")
   return(paste("D2Psi.d.",  - J, ".", filter.number, ".", family,sep = ""))
if(switch == "level")
   return(paste("D2Psi.l.",  - J, ".", filter.number, ".", family,sep = ""))
else 
   stop("\nSwitch can only take the values direction and level!\n")
}

