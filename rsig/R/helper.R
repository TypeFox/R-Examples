filterSD = function(X, sd.filter) {
  if (exists("topk", sd.filter)) {
    return(rank(apply(X, 2, sd), ties.method = "random") <= sd.filter$topk)
  } else if (exists("quant", sd.filter)) {
    tmp = apply(X, 2L, sd)
    return(tmp >= quantile(tmp, sd.filter$quant))
  } 
  stop("filterSD: Either topk or quant must be passed in sd.filter")
}
