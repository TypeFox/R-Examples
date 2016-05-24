# This file contains functions to call phyclust_find_consensus() in C.

find.consensus <- function(X, code.type = .code.type[1], with.gap = FALSE){
  code.type <- which(code.type == .code.type) - 1

  ret <- .Call("R_phyclust_find_consensus",
               as.integer(nrow(X)),
               as.integer(ncol(X)),
               as.integer(code.type),
               as.integer(with.gap),
               as.integer(t(X)),
               PACKAGE = "phyclust")
  ret
} # End of find.consensus().

