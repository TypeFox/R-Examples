#' Creates the negative log-likelihood function
#' \code{create_nll} Creates the negative log-likelihood function
#' @keywords internal
#' @export
create_nll <- function(d, x, k, n, psyfunguesslapses){
  function(p) {
    phi <- psyfunguesslapses(d[[x]], p)
    phi[phi < .Machine$double.eps] <- .Machine$double.eps
    phi[phi > (1 - .Machine$double.eps)] <- 1 - .Machine$double.eps
    return(-sum(d[[k]] * log(phi) + (d[[n]] - d[[k]]) * log(1 - phi)))
  }
}


