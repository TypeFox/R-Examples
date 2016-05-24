#' @rdname full2meta
#' @param d The dimension of the mixture distribution.
#' @return meta2full returns a formatted list of parameters as described by
#'   \code{\link{rtheta}}.
#' @export
meta2full <- function(par, d) {  # par = c(pie1, mu, sigma, rho)
  if (par[4] <= -1/(d-1))
    stop("correlation coefficient is not valid")
  list(m = 2, d = d, pie = c(pie1 = unname(par[1]), pie2 = unname(1-par[1])),
       mu = list(comp1 = rep(0, d), comp2 = rep(par[2], d)),
       sigma = list(comp1 = diag(d),
                    comp2 = matrix(par[4]*par[3]^2, d, d) +
                      diag(rep((1-par[4])*par[3]^2, d))))
}

