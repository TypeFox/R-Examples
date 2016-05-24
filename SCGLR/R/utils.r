# weighted scale functions
# @param x vector to scale
# @param w weight
# @return scaled vector
wtScale <-function(x,w) {
  xc=x-sum(w*x)
  v=sum(xc^2*w)
  xcr=xc/sqrt(v)
  return(xcr)                                   
}

# utils functions
# @param x vector to center
# @param w weight
# @return centered vector
wtCenter=function(x,w) {
  xc=x-sum(w*x)
  return(xc)
}

checkLossFunction <- function(type) {
  if(!type %in% c("auc","likelihood","aic","aicc","bic","mspe"))
    stop("Unknown loss function!")
}

# returns string w/o leading or trailing whitespace
trim <- function (x) gsub("^\\s+|\\s+$", "", x)
