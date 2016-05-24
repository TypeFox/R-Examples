pglpenalty <-
function(X,p,b,type = c("growth","decay","both"),plotpenalty = TRUE, allowed.error = 0.005,invert = FALSE) {
  type = match.arg(type)
  
  if (b<=0) {
    stop("b (rate of growth) must be greater than 0")
  } else if (p <= 0 || p>= 1) {
    stop("p,percent of vector subject to penalty, must be greater than 0 and less than 1")
  } else if (p > 0.5 & type == "both") {
    stop("p cannot be greater than 0.5 if type = both")
  }
  
  i1 = as.numeric(X[ceiling(length(X)*p)])
  M = (log((1/(1 - allowed.error))-1)/b) + i1
  i0 = (log((1/allowed.error) - 1)/-b) + M
  
  if (type == "both") {
    x1 = c(i1,max(X) - i1)
    x0 = c(i0,x1[2] + (i1-i0))
  } else if (type == "growth") {
    x1 = i1
    x0 = i0
  } else if (type == "decay") {
    x1 = max(X) - i1
    x0 = x1 + (i1-i0)
  }
  pglpenalty = vglpenalty(X,x0,x1,plotpenalty,allowed.error,invert)
  pglpenalty
}
