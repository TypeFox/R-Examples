print.aodql <- function(x, ...){
  res <- list(model = x$fm,
              method = x$method,
              phi = x$phi,
              phi.scale = x$phi.scale,
              Z = x$Z,
              nbiter = x$nbiter)
  print(res, ...)
  invisible(res)
  }