#Internal Function
acts2probs <- function(acts) {
  acts.minimum.correction=1e-10
  acts.minimum = abs(min(acts))
  if(acts.minimum==0)
    acts.minimum = acts.minimum.correction
  acts = acts + acts.minimum

  rowsums = apply(acts,1,sum)

  m = matrix(rep(rowsums, rep(ncol(acts),nrow(acts))), nrow(acts), ncol(acts), byrow=TRUE)

  p = acts/m
  colnames(p) = colnames(acts)

  predictions = apply(p, 1, function(v) names(which.max(v)))

  result <- list(p = p, predicted = predictions)
  return(result)
}

