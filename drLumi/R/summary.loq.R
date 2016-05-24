#' @export
summary.loq <- function(object, ... ) {
  x <- object
  if (!inherits(x,"loq")) {
    stop("not loq object")
  }
  
  lloq <- unlist(lapply(x, function(y) y$lloq))
  uloq <- unlist(lapply(x, function(y) y$uloq))
  ly <- unlist(lapply(x, function(y) y$ly))
  uy <- unlist(lapply(x, function(y) y$uy))
  loq.drange <- abs(uloq - lloq)
  y.drange <- abs(uy - ly)
  method <- unlist(lapply(x, function(y) y$method))
  
  ans <- data.frame(analyte = names(x), lloq=lloq, uloq=uloq, ly=ly, uy=uy,
                    loq.drange = loq.drange,
                    y.drange = y.drange,
                    method=method)
  rownames(ans) <- 1:nrow(ans)
  return(ans)
}
