#' @rdname SimulateGMCMData
#' @return \code{SimulateGMMData} returns a list of length 3 with elements:
#'   \item{z}{A matrix of GMM realizations.}
#'   \item{K}{An integer vector denoting the component from which the
#'     realization comes.}
#'   \item{theta}{As above and in \code{\link{rtheta}}.}
#' @export
SimulateGMMData <- function(n = 1000, theta = rtheta(...), ...) {

  K <- sample(seq_len(theta$m), size = n, replace = TRUE, prob = theta$pie)

  if (length(unique(K)) != theta$m) {
    warning(paste("Some components were not represented in the sample.",
                  "The mixture proportions may be too small to be sampled",
                  "with the given n.",
                  "Ignore this warning, try to resample, or change theta."))
  }

  tab <- table(K)

  SampFunc <- function (i) {
    comp.k <- as.numeric(names(tab[i]))

    # Using rmvnormal instead of mvtnorm::rmvnorm
    rmvnormal(tab[i],
              mu = theta$mu[[comp.k]],
              sigma = theta$sigma[[comp.k]])

    ## Not using mvtnorm anymore.
    #mvtnorm::rmvnorm(tab[i], theta$mu[[comp.k]], theta$sigma[[comp.k]])
  }

  z <- do.call(rbind, lapply(seq_along(tab), FUN = SampFunc))[order(order(K)), ]
  return(list(z = unname(z), K = K, theta = theta))
}
