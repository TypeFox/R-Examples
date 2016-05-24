predict.CCM <-
function(object, y, func = mean, ret.scores = FALSE, ...) {
  keep = !is.na(y)
  y = y[keep]
  object = object[,keep]
  preds = rep(NA,nrow(object)) 
  scores = NULL
  for (i in 1:nrow(object)) {
    s = split(object[i,], y)
    l = (lapply(s,func, na.rm=TRUE, ...))
    w = which.max(l)
    if (ret.scores) scores = rbind(scores, unlist(l))
    preds[i] = names(w)
  }
  if (ret.scores) {
        scores = t(scores)
	return(scores)
  }
  if (is.numeric(y)) preds = as.numeric(preds)
  if (is.factor(y)) preds = as.factor(preds)

  return(preds)
}

