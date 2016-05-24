list2matrix <- function(x, diag=FALSE) {
  if (!is.list(x))
    stop("\"x\" has to be a list.")
  if (!identical(0, var(sapply(x, function(x){dim(x)[[1]]}))))
    stop("Dimensions of matrices in \"x\" have to be the same in order to stack them together.")
  
  if (is.null(dimnames(x[[1]]))) {
    oldNames <- paste("x", 1:dim(x[[1]])[[1]], sep = "")
  } else {
    oldNames <- dimnames(x[[1]])[[1]]
  }
   
  if (diag) {
    psNames <- vech(outer(oldNames, oldNames, paste, sep = "_"))
    out <- t(sapply(x, function(x) {(vech(x))}))
  } else {
    psNames <- vechs(outer(oldNames, oldNames, paste, sep = "_"))
    out <- t(sapply(x, function(x) {(vechs(x))}))
  }

  dimnames(out) <- list(names(x), psNames)
  out
}
