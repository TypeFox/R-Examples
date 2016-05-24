#' @title Compute (weighted) covariance matrix from frequency spectrum

#' @description
#' \code{mvspectrum2wcov} computes a (weighted) covariance matrix estimate 
#' from the frequency spectrum (see Details).
#' 
#' @details
#' The covariance matrix of a multivariate time series satisfies the identity
#' \deqn{
#' \Sigma_{X} \equiv \int_{-\pi}^{\pi} S_{X}(\lambda) d \lambda.
#' }
#' 
#' A generalized covariance matrix estimate can thus be obtained using a weighted average
#' \deqn{
#' \tilde{\Sigma}_X = \int_{-\pi}^{\pi} K(\lambda) S_{X}(\lambda) d \lambda,
#' }
#' where \eqn{K(\lambda)} is a kernel symmetric around \eqn{0} which averages out to
#' \eqn{1} over the interval \eqn{[-\pi, \pi]}, i.e.,
#'  \eqn{\frac{1}{2 \pi} \int_{-\pi}^{\pi} K(\lambda) d \lambda = 1}. 
#'  This allows one to remove or amplify specific frequencies in the covariance matrix 
#'  estimation.
#' 
#' For ForeCA \code{mvspectrum2wcov} is especially important as we use
#' \deqn{
#' K(\lambda) = -\log f_y(\lambda),
#' }
#' as the \emph{weights} (their average is not \eqn{1}!).  This particular kernel
#' weight is implemented as a wrapper in \code{weightvector2entropy_wcov}.
#' 
#' @inheritParams common-arguments
#' @param kernel.weights numeric; weights for each frequency. By default uses 
#' weights that average out to \code{1}.
#' @return 
#' A symmetric \eqn{n \times n} matrix.
#' 
#' If \code{kernel.weights} \eqn{\geq 0}, then it is positive semi-definite;
#' otherwise, it is symmetric but not necessarily positive semi-definite.
#' @seealso \code{\link{mvspectrum}}
#' @keywords ts
#' @export
#' @examples
#' 
#' nn <- 50
#' YY <- cbind(rnorm(nn), arima.sim(n = nn, list(ar = 0.9)), rnorm(nn))
#' XX <- YY %*% matrix(rnorm(9), ncol = 3)  # random mix
#' XX <- scale(XX, scale = FALSE, center = TRUE)
#' 
#' # sample estimate of covariance matrix
#' Sigma.hat <- cov(XX)
#' dimnames(Sigma.hat) <- NULL
#' 
#' # using the frequency spectrum
#' SS <- mvspectrum(XX, "wosa")
#' Sigma.hat.freq <- mvspectrum2wcov(SS)
#' 
#' layout(matrix(1:4, ncol = 2))
#' par(mar = c(2, 2, 1, 1))
#' plot(c(Sigma.hat/Sigma.hat.freq))
#' abline(h = 1)
#' 
#' image(Sigma.hat)
#' image(Sigma.hat.freq)
#' image(Sigma.hat / Sigma.hat.freq)
#' 

mvspectrum2wcov <- function(mvspectrum.output, kernel.weights = 1) {
  nseries <- 1
  if (!is.null(dim(mvspectrum.output))) {
    nseries <- dim(mvspectrum.output)[2]
  }
  
  if (nseries == 1) {
    mvspectrum.output <- array(mvspectrum.output, 
                               c(length(mvspectrum.output), 1, 1))
  }
  
  # change the function to be integrated if kernel.weights != 1
  if (any(kernel.weights != 1)) {
    mvspectrum.output <- sweep(mvspectrum.output, 1, kernel.weights, "*")
  } 
  weighted.cov <- apply(mvspectrum.output, 2:3, sum)
  weighted.cov <- 2 * Re(weighted.cov)  # same as: weighted.cov + Conj(weighted.cov)
  return(weighted.cov)
} 

#' @rdname mvspectrum2wcov
#' @export
#' @description
#' \code{weightvector2entropy_wcov} computes the weighted covariance
#' matrix using the negative entropy of the univariate spectrum (given the
#' weightvector) as kernel weights.  This matrix is the objective matrix
#' for many \code{foreca.*} algorithms.
#' @inheritParams foreca.EM.h
#' @inheritParams foreca.EM.E_step
#' @examples
#' # examples for entropy wcov
#' XX <- diff(log(EuStockMarkets)) * 100
#' UU <- whiten(XX)$U
#' ff <- mvspectrum(UU, 'wosa', normalize = TRUE)
#' 
#' ww0 <- initialize_weightvector(num.series = ncol(XX), method = 'rnorm')
#' 
#' weightvector2entropy_wcov(ww0, ff,
#'                           entropy.control = 
#'                             list(prior.weight = 0.1))

weightvector2entropy_wcov <- function(weightvector = NULL, f.U, 
                                      f.current = NULL,
                                      entropy.control = list()) {
    
    check_mvspectrum_normalized(f.U)
    num.freqs <- dim(f.U)[1]
    num.series <- dim(f.U)[2]
    
    entropy.control <- complete_entropy_control(entropy.control,
                                                num.outcomes = 2 * num.freqs)
    
    if (entropy.control$prior.weight > 0) {
      mix.weight <- entropy.control$prior.weight
      # divide by length so that when summing over all frequencies we
      # get original weight back
      mix.weight <- mix.weight / num.freqs
      
      # make a prior 3D spectrum with frequencies in first and diagonal
      # matrix with prior weight in diagonal matrix of 2nd and 3rd dimension
      prior.3d <- entropy.control$prior.probs[seq_len(num.freqs)]
      dim(prior.3d) <- c(num.freqs, 1, 1)
      prior.3d <- apply(prior.3d, 1, function(x) diag(as.numeric(x), num.series))
      dim(prior.3d) <- c(num.series, num.series, num.freqs)
      prior.3d <- base::aperm(prior.3d, perm = c(3, 1, 2))
      #prior.3d <- R.utils::wrap(prior.3d, map = list(3, 1, 2))
      
      # check_mvspectrum_normalized(prior.3d)
      # f.current.prior <- spectrum_of_linear_combination(prior.3d,
      #                                                           weightvector)
      # check_mvspectrum_normalized(f.current.prior)
      stopifnot(dim(prior.3d) == dim(f.U))
      f.U <- (1 - mix.weight) * f.U + mix.weight * prior.3d
    }
    
    if (is.null(f.current)) {
      f.current <- spectrum_of_linear_combination(f.U, weightvector)
    }
    check_mvspectrum_normalized(f.current)
    
    # keep only positive values since 0 * log(0) = 0 in entropy calculation
    f.U <- f.U[f.current > 0,,]
    f.current <- f.current[f.current > 0]
    
    # don't check f.U again since the smoothing does not guarantee 0 off diagonal
    check_mvspectrum_normalized(f.current)
    # dont normalize f.U necessarily, since f.current = 0, not f.U
    # multiply by num.freqs for numerical stability
    SS <- mvspectrum2wcov(num.freqs * f.U, kernel.weights = 
                            -log(f.current, base = entropy.control$base))
    # divide by num.freqs again
    SS <- SS / num.freqs
    return(SS)
  }
