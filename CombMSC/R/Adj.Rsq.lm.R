`Adj.Rsq.lm` <-
function(object, ...){
  -1 * summary(object)$adj.r.squared
}

