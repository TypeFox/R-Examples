genComp = function(n, len = TRUE, addZeros = FALSE){
  regularize = function(x, regLen){
    return(c(x, rep(0, regLen - length(x))))
  }
  if(len & is.logical(len)){
    
    l = generateCompositions(n)
    if(addZeros)
      l = lapply(l, regularize, regLen = n)
    return(l)
    
  }else if(is.numeric(len) & len <= n) {
    
    constraintFn = function(l){
      if(length(l) <= len)
        return(TRUE)
      else
        return(FALSE)
    }
    
    l = generateCompositions(n)
    l = l[sapply(l, constraintFn)]
    
    if(addZeros)
      l = lapply(l, regularize, regLen = len)
    
    return(l)
  }
}