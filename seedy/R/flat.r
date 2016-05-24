flat <- function(x, span, eq.size=5000) {
  if (eq.size%%1!=0 || eq.size<=0) {
    stop("Equilibrium population size must be a postive integer")
  }
  return(eq.size)
}