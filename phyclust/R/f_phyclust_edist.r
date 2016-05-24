### This file contains functions for compute evolution distances.

phyclust.edist <- function(X, edist.model = .edist.model[1]){
  if(is.vector(X)) stop("The X should be a matrix and nrow > 1.")

  if(edist.model[1] %in% .edist.model){
    edist.model <- which(edist.model[1] == .edist.model) - 1
  } else{
    stop("The distance model is not found.")
  }
  N <- nrow(X)
  L <- ncol(X)

  d <- .Call("R_phyclust_edist",
             as.integer(edist.model),
             as.integer(N),
             as.integer(L),
             as.integer(t(X)),
             PACKAGE = "phyclust")
  class(d) <- "dist"
  attr(d, "Size") <- N
  attr(d, "Diag") <- FALSE
  attr(d, "Upper") <- FALSE
  attr(d, "method") <- .edist.model[edist.model[1] + 1]
  d
} # End of edist().

