`A2name` <-
function(J, filter.number, family, switch = "direction")
{
#
# Program to return a specific character string format
# for a given 2-D inner product matrix of 
# discrete autocorrelation wavelets
#
if(J >= 0.) stop("J must be a negative integer")
if(switch == "direction") {
   return(paste("D2Amat.d.",  - J, ".", filter.number, ".", family, sep = ""))
}
if(switch == "level") {
   return(paste("D2Amat.l.",  - J, ".", filter.number, ".", family, sep = ""))
}
else {
  stop("You can only allow switch to equal direction or level!"
)
}
}

