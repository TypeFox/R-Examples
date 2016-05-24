######################################################################
# Take any confidence interval generating function for the binomial
# proportion and increase n until a solution is found
# which has the desired width
######################################################################

ciss.binom <- function(p0, d, alpha=0.05, ci.fun=binom.confint, np02x = function(n, p0) round(n*p0), verbose=FALSE,nStart=1,nMax=1e6,...) {
  if (nStart < 0) {
    stop("nStart has to be an integer >= 0.")
  }
  
  n <- max(0,round(nStart - 1))
  stop <- FALSE

  #Loop until a solution is found where the interval width is below 2*d
  while (!stop & (n<nMax)) {
    n <- n+1

    #Check if ... contains arguments to pass to np02x function
    args <- list(n=n, p0=p0)
    extras <- match.call(expand.dots=FALSE)$...
    if (length(extras) > 0) {
       existing <- !is.na(match(names(extras), names(formals(np02x))))
       for (a in names(extras)[existing]) {
         args[[a]] <- extras[[a]]
       }
    }
    #Find corresponding x value leading to pi0 estimate by calling np02x    
    x <- do.call("np02x",args )
    #Compute width of corresponding confidence interval
    ci <- ci.fun(x=x, n=n, conf.level=1-alpha,...)[,c("lower","upper")]
    width <- as.numeric(ci[2]-ci[1])
    if (verbose) {
      cat("n = ",n," width = ", width, "\n")
    }
    stop <- width < 2*d
  }
  if (n==nMax) {
    stop("No solution for n found which is smaller than nMax")
  }
  return(n)
}
  
