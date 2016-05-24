coef.mvdareg <- function (object, ncomp = object$ncomp, 
  type = c("coefficients", "loadings", "weights", "y.loadings"), conf = .95, ...) {
  switch(match.arg(type), weights = {
    weights.mvdareg(object, ncomp, conf, ...)
  }, 
  loadings = {
    loadings.mvdareg(object, ncomp, conf, ...)
  }, 
  coefficients = {
    coefficients.mvdareg(object, ncomp, conf, ...)
  }, 
  y.loadings = {
    y.loadings(object, conf, ...)
  })
}
