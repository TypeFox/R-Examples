bglpenalty <-
function(X,x1,b,type = c("growth","decay","both"),plotpenalty = TRUE, allowed.error = 0.005,invert = FALSE) {
  type = match.arg(type)
  
  if (b<=0) {
    stop("b (rate of growth) must be greater than 0")
  } else if (length(x1)> 2) {
    stop("vector x1 must be contain only 1 or 2 numbers")
  }
  
  if (length(x1)==1) {
    i1 = x1
    M = (log((1/(1 - allowed.error))-1)/b) + i1
    i0 = (log((1/allowed.error) - 1)/-b) + M
    
    if (type == "both") {
      x1 = c(i1,max(X) - i1)
      x0 = c(i0,x1[2] + (i1-i0))
    } else if (type == "growth") {
      x1 = i1
      x0 = i0
    } else if (type == "decay") {

      x0 = x1 + (i1-i0)
    }
    
  } else {
    i1 = x1[1]
    M = (log((1/(1 - allowed.error))-1)/b) + i1
    i0 = (log((1/allowed.error) - 1)/-b) + M
    x0 = c(i0,x1[2] + (i1-i0))
  }
  
  gglpenalty = vglpenalty(X,x0,x1,plotpenalty,allowed.error,invert)
  gglpenalty
}
