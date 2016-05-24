hump <- function(x, span, max.size=5000) {
  if (max.size%%1!=0 || max.size<=0) {
    stop("Equilibrium population size must be a postive integer")
  }
  return(max.size*(cos(2*pi*x/span-pi)+1)/2)
}