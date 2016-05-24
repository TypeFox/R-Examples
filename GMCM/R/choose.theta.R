#' Heuristically chosen starting value of theta
#'
#' This function uses a \code{k}-means algorithm to heuristically select
#' suitable starting values for the general model.
#'
#' The function selects the centers from the k-means algorithm as an initial
#' estimate of the means. The proportional sizes of the clusters are selected
#' as the initial values of the mixture proportions. The within cluster
#' standard deviations are used as the variance of the clusters. The
#' correlations between each dimension are taken to be zero.
#'
#' @param u A matrix of (estimates of) realizations from the GMCM.
#' @param m The number of components to be fitted.
# @param fac Numeric. A factor applied to the standard deviation estimates of
#   each component and each dimension.
#' @param \dots Arguments passed to \code{\link{kmeans}}.
#' @return A list of parameters for the GMCM model on the form described in
#'   \code{\link{rtheta}}.
#' @note The function uses the \code{kmeans} function from the
#'   \code{stats}-package.
#' @author Anders Ellern Bilgrau <anders.ellern.bilgrau@@gmail.com>
#' @examples
#' set.seed(2)
#'
#' # Simulating data
#' data1 <- SimulateGMCMData(n = 10000, m = 3, d = 2)
#' obs.data <- Uhat(data1$u)  # The ranked observed data
#'
#' # Using choose.theta to get starting estimates
#' theta <- choose.theta(u = obs.data, m = 3)
#' print(theta)
#'
#' # To illustrate theta, we simulate from the model
#' data2 <- SimulateGMMData(n = 10000, theta = theta)
#'
#' cols <- apply(get.prob(obs.data,theta),1,which.max)
#'
#' # Plotting
#' par(mfrow = c(1,3))
#' plot(data1$z, main = "True latent GMM")
#' plot(Uhat(data1$u), col = cols,
#'      main = "Observed GMCM\nColoured by k-means clustering")
#' plot(data2$z, main = "initial GMM")
#' @export
choose.theta <- function(u,
                         m,
                         #fac = 2,
                         ...) {

  km  <- kmeans(u, centers = m, ...)
  pie <- km$size/sum(km$size)  #  Estimate of pie
  mu  <- lapply(1:m, function(i) km$centers[i, ])

  get.sigma <- function(i) {
    diag(colSds(u[km$cluster == i, , drop = FALSE]))/sqrt(km$size[i])
  }
  sigma <- lapply(1:m, get.sigma)
  #sigma <- lapply(sigma, "*", fac^2)

  # Scaling and translating
  mu      <- lapply(mu, "-", mu[[1]])           # Translating means

  scaling <- diag(sigma[[1]])[1]                # Get scaling factor
  sigma   <- lapply(sigma, "/", scaling)        # Scaling sigma
  mu      <- lapply(mu, "/", sqrt(scaling))     # Scaling means

  # Naming
  names(pie) <- paste("pie", 1:m, sep = "")
  names(mu) <- names(sigma) <- paste("comp", 1:m, sep = "")
  ans <- list(m = m, d = ncol(u), pie = pie, mu = mu, sigma = sigma)

  return(ans)
}
