cvpvs <-
function(X, Y, method = c('gaussian', 'knn', 'wnn', 'logreg'), ...){
  method <- match.arg(method)
  pv <- switch(method,
               gaussian = cvpvs.gaussian(X = X, Y = Y, ...),
               knn = cvpvs.knn(X = X, Y = Y, ...),
               wnn = cvpvs.wnn(X = X, Y = Y, ...),
               logreg = cvpvs.logreg(X = X, Y = Y, ...))
  return(pv)
}
