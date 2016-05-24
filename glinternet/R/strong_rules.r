strong_rules = function(candidates, lambda, lambdaPrev){

  constant = 2*lambda - lambdaPrev
  if (constant <= 0) return(candidates$variables)
  
  labels = c("cat", "cont", "catcat", "contcont", "catcont")
  
  strongSet = lapply(1:5, function(x) if(!is.null(candidates$variables[[labels[x]]])){
    index = candidates$norms[[labels[x]]] >= constant
    if (any(index)) matrix(candidates$variables[[labels[x]]][index, ], nrow=sum(index))
    else NULL
  })
  names(strongSet) = labels

  return (strongSet)
}
