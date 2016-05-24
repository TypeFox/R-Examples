print.multi_tpfit <-
function(x, ...) {
  nc <- length(x$coefficients)
  dire.mat <- diag(, nc)
  if (!is.null(x$rotation)) {
    dire.mat <- .C('rotaxes', nc = as.integer(nc), ang = as.double(x$rotation),
                   res = as.double(dire.mat), PACKAGE = "spMC")$res
    dire.mat <- t(matrix(dire.mat, nc, nc))
  } 
  for(i in 1:nc) {
    cat("Direction (", sep = "")
    cat(dire.mat[i,], sep = ", ")
    cat(")\n", sep = "")
    print(x$coefficients[[i]], ...)
  }
  cat("Estimated proportions:\n\n")
  proportions <- as.vector(x$prop)
  names(proportions) <- names(x$prop)
  print(proportions, ...)
  cat("\n")
  invisible(x)
}

