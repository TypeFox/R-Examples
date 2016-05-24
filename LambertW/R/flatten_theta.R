#' @rdname theta-utils
#' @description
#' \code{flatten_theta} and \code{unflatten_theta} convert between the list 
#' \code{theta} and its vector-style flattened type.  The flattened version is required
#' for several optimization routines, since they optimize over multivariate vectors -- not lists.
#' @export
flatten_theta <- function(theta) {
  
  params <- unlist(theta)
  # if theta was emtpy, then params = NULL; no name changes are necessary (would result in error)
  if (!is.null(params)) {
    # remove 'beta.' prefix
    names(params) <- gsub("^beta.", "", names(params))
    # remove 'alpha.', 'gamma.', 'delta.' prefix
    for (var.name in c("alpha", "gamma", "delta")) {
      names(params) <- gsub(paste0("^", var.name, "."), "", names(params))
    }
  }
  return(params)
} 
