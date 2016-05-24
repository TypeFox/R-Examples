pvs <-
function(NewX, X, Y, method = c('gaussian', 'knn', 'wnn', 'logreg'), ...) {
  method <- match.arg(method)
  pv <- switch(method,
               gaussian = pvs.gaussian(NewX = NewX, X = X, Y = Y, ...),
               knn = pvs.knn(NewX = NewX, X = X, Y = Y, ...),
               wnn = pvs.wnn(NewX = NewX, X = X, Y = Y, ...),
               logreg = pvs.logreg(NewX = NewX, X = X, Y = Y, ...))
  return(pv)
}
