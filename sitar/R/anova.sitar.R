anova.sitar <- function (object, ..., test = TRUE, type = c("sequential", "marginal"),
                         adjustSigma = TRUE, Terms, L, verbose = FALSE) {
  name <- deparse(substitute(object))
  class(object) <- class(object)[class(object) != 'sitar']
  assign(name, object)
  dots <- match.call(expand.dots=FALSE)$...
  if (length(dots) > 0) for (obj in dots) {
    name <- deparse(obj)
    obj <- eval(obj)
    class(obj) <- class(obj)[class(obj) != 'sitar']
    assign(name, obj)
  }
  do.call('anova', as.list(match.call()[-1]))
}
