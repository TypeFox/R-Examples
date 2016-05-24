#' @method clean.scales cdslist
#' @rdname clean.scales
#' @rdname clean.scales
#' @export
clean.scales.cdslist <- function(object, data, K, col.subset = NULL, ...){
  Kvec <- lapply(object, "[[", "K")
  object <- object[[which(Kvec == K)]]
  return(clean.scales(object, data, col.subset, ...))
}