#' @rdname theta-utils
#' @param theta.flattened named vector; flattened version of list \code{theta}.
#' @export
unflatten_theta <- function(theta.flattened, distname, type) {

  check_distname(distname)

  # where are parameters that start with alpha, gamma, or delta
  not.beta.pos <- grepl("^alpha|^gamma|^delta", names(theta.flattened))
  theta <- list(beta = theta.flattened[!not.beta.pos])
  
  if (type == "s") {
    theta[["gamma"]] <- theta.flattened["gamma"]
    names(theta[["gamma"]]) <- NULL
  } else {
    alpha.ind <- grepl("^alpha", names(theta.flattened))
    if (any(alpha.ind)) {
      theta[["alpha"]] <- theta.flattened[alpha.ind]
    }
    delta.ind <- grepl("^delta", names(theta.flattened))
    if (any(delta.ind)) {
      theta[["delta"]] <- theta.flattened[delta.ind]
    }
  }
  return(theta)
}
