summary.diagram <-
function(object, ...) {
  call <- attributes(object)[["call"]]
  n <- nrow(object)
  maxdimension <- attributes(object)[["maxdimension"]]
  scale <- attributes(object)[["scale"]]
  out <- list(call = call, n = n, maxdimension = maxdimension, scale = scale)
  class(out) <- c("summary.diagram")
  return(out)
}