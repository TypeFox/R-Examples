predict.OKriglistplus <- function(object, X, ...) {
  ## serves eg in sequence safeSurface->surface.OKrig->predict_surface->predict
  ## if the object is a OKriglistplus (i.e. predict is this predict.Kriglistplus).
  ## note that predict.OKriglistplus->purefn->predict.Krig.
  ## X can be a matrix (vector of parameter points).
  nrX <- nrow(X)
  if (is.null(nrX)) {
    llocalst <- "(!) From predict.Kriglistplus(..., X, ...): nrow(X) is NULL. For single points, use purefn()."
  }
  rv <- apply(X, 1, purefn, fitobj=object, ...)
  return(rv)
}
