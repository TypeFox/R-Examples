##' @export
parlabels <- function(x,exo=FALSE) {
  res <- c(unlist(intfix(x)[unlist(lapply(intfix(x), function(y) !is.na(y) & !is.numeric(y)))]),
           regfix(x)$labels[!is.na(regfix(x)$labels)],
           covfix(x)$labels[!is.na(covfix(x)$labels)])
  if (!is.null(x$exfix))
    res <- c(res,
           unlist(x$exfix[!is.na(x$exfix) && !is.numeric(x$exfix)]))
  if (exo)
    res <- intersect(res,index(Model(x))$exogenous)
  return(res)
}
