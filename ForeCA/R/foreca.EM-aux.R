#' @title ForeCA EM auxiliary functions
#' @name foreca.EM-aux
#' @description 
#' \code{foreca.EM.one_weightvector} relies on several auxiliary functions:
#' 
NULL

#' @rdname foreca.EM-aux
#' @description
#' \code{foreca.EM.E_step} computes the spectral density of 
#' \eqn{y_t=\mathbf{U}_t \mathbf{w}} given the weightvector \eqn{\mathbf{w}} 
#' and the normalized spectrum estimate \eqn{f_{\mathbf{U}}}.
#' A wrapper around \code{\link{spectrum_of_linear_combination}}.
#' @keywords manip math
#' 
#' @inheritParams common-arguments
#' @param weightvector numeric; weights \eqn{\mathbf{w}} for 
#' \eqn{y_t = \mathbf{U}_t \mathbf{w}}. Must have unit norm in \eqn{\ell^2}. 
#' @export
#' @return
#' \code{foreca.EM.E_step} returns the normalized univariate spectral 
#' density (normalized such that its \code{sum} equals \eqn{0.5}).
#' @examples
#' XX <- diff(log(EuStockMarkets)) * 100
#' UU <- whiten(XX)$U
#' ff <- mvspectrum(UU, 'wosa', normalize = TRUE)
#' 
#' ww0 <- initialize_weightvector(num.series = ncol(XX), method = 'rnorm')
#' 
#' f.ww0 <- foreca.EM.E_step(ff, ww0)
#' plot(f.ww0, type = "l")

foreca.EM.E_step <- function(f.U, weightvector) {
  
  stopifnot(!any(is.na(weightvector)),
            isTRUE(all.equal(target = 1, 
                             current = base::norm(weightvector, "2"))),
            length(dim(f.U)) == 3,
            dim(f.U)[2] == dim(f.U)[3], # square
            dim(f.U)[2] == length(weightvector))  # weightvector has the right dimension
  
  check_mvspectrum_normalized(f.U)
  
  spec.dens.est <- spectrum_of_linear_combination(f.U, weightvector)
  check_mvspectrum_normalized(spec.dens.est)  # is also a test already
  
  return(spec.dens.est)
} 

#' @rdname foreca.EM-aux
#' @description
#' \code{foreca.EM.M_step} computes the minimizing eigenvector 
#' (\eqn{\rightarrow \widehat{\mathbf{w}}_{i+1}}) of the weighted
#' covariance matrix, where the weights equal the negative logarithm of the 
#' spectral density at the current \eqn{\widehat{\mathbf{w}}_i}.
#' @keywords manip
#' @inheritParams common-arguments
#' @param minimize logical; if \code{TRUE} (default) it returns the eigenvector
#' corresponding to the smallest eigenvalue; otherwise to the largest eigenvalue.
#' @inheritParams common-arguments
#' @export
#' @return
#' \code{foreca.EM.M_step} returns a list with three elements:
#' \itemize{
#'    \item \code{matrix}: weighted covariance matrix, where the weights are 
#'          the negative log of the spectral density.  If density is estimated 
#'          by discrete probabilities, 
#'          then this \code{matrix} is positive semi-definite, since 
#'          \eqn{-\log(p) \geq 0} for \eqn{p \in [0, 1]}. 
#'          See \code{\link{weightvector2entropy_wcov}}.
#'    \item \code{vector}: minimizing (or maximizing if 
#'          \code{minimize = FALSE}) eigenvector of \code{matrix},
#'    \item \code{value}: corresponding eigenvalue.
#'    }
#' @examples
#' 
#' one.step <- foreca.EM.M_step(ff, f.ww0, 
#'                              entropy.control = list(prior.weight = 0.1))
#' image(one.step$matrix)
#' \dontrun{
#' requireNamespace(LICORS)
#' # if you have the 'LICORS' package use
#' LICORS::image2(one.step$matrix)
#' }
#' ww1 <- one.step$vector
#' f.ww1 <- foreca.EM.E_step(ff, ww1)
#' 
#' layout(matrix(1:2, ncol = 2))
#' matplot(seq(0, pi, length = length(f.ww0)), cbind(f.ww0, f.ww1), 
#'         type = "l", lwd =2, xlab = "omega_j", ylab = "f(omega_j)")
#' plot(f.ww0, f.ww1, pch = ".", cex = 3, xlab = "iteration 0", 
#'      ylab = "iteration 1", main = "Spectral density")
#' abline(0, 1, col = 'blue', lty = 2, lwd = 2)
#' 
#' Omega(mvspectrum.output = f.ww0) # start
#' Omega(mvspectrum.output = f.ww1) # improved after one iteration

foreca.EM.M_step <- function(f.U, f.current, minimize = TRUE,
                             entropy.control = list()) {
  stopifnot(all(f.current >= 0),
            is.logical(minimize))
  
  num.freqs <- length(f.current)
  
  check_mvspectrum_normalized(f.U)
  check_mvspectrum_normalized(f.current)
  
  entropy.control <- complete_entropy_control(entropy.control,
                                              num.outcomes = 2 * num.freqs)
  integrated.spectrum.entropy <- 
      weightvector2entropy_wcov(NULL, 
                                f.U = f.U,
                                f.current = f.current,
                                entropy.control = entropy.control)
  EE <- eigen(integrated.spectrum.entropy, symmetric = TRUE)
  sel <- ifelse(minimize, which.min(EE$values), which.max(EE$values))
  
  weightvector <- EE$vector[, sel]
  
  # make first entry always positive (scale vector always this way)
  # for consistent results among re-runs
  weightvector <- weightvector / sign(weightvector[1])
  out <- list(matrix = integrated.spectrum.entropy, 
              vector = weightvector, 
              value = EE$values[sel])
  return(out)
} 

#' @rdname foreca.EM-aux
#' @description
#' \code{foreca.EM.E_and_M_step} is a wrapper around \code{foreca.EM.E_step}
#' followed by \code{foreca.EM.M_step}.
#' @keywords manip
#' @inheritParams common-arguments
#' @export
#' @return
#' Contrary to \code{foreca.EM.M_step}, \code{foreca.EM.E_and_M_step} only returns the optimal 
#' weightvector as a numeric. 
#' @examples
#' 
#' ww0 <- initialize_weightvector(NULL, ff, method = "rnorm")
#' ww1 <- foreca.EM.E_and_M_step(ww0, ff)
#' ww0
#' ww1
#' barplot(rbind(ww0, ww1), beside = TRUE)
#' abline(h = 0, col = "blue", lty = 2)
#' 

foreca.EM.E_and_M_step <- function(weightvector, f.U, minimize = TRUE,
                                   entropy.control = list()) {
  
  f.current <- spectrum_of_linear_combination(f.U, weightvector)
  opt.vec <- foreca.EM.M_step(f.U, f.current, minimize = minimize,
                              entropy.control = entropy.control)$vector
  return(opt.vec)
} 


#' @rdname foreca.EM-aux
#' @description
#' \code{foreca.EM.h} evaluates (an upper bound of) the entropy of the spectral density as a function
#' of \eqn{\mathbf{w}_i} (or \eqn{\mathbf{w}_{i+1}}). This is the objective funcion that should be 
#' \code{minimize}d.
#' 
#' @keywords manip
#' @param weightvector.new weightvector \eqn{\widehat{\mathbf{w}}_{i+1}} of the new 
#' iteration (i+1).
#' @param f.current numeric; spectral density estimate of 
#' \eqn{y_t=\mathbf{U}_t \mathbf{w}} for the current estimate 
#' \eqn{\widehat{\mathbf{w}}_i} (required for 
#' \code{foreca.EM.M_step}; optional for \code{foreca.EM.h}).
#' @param weightvector.current weightvector \eqn{\widehat{\mathbf{w}}_{i}} of the 
#' current iteration (i).
#' @param return.negative logical; if \code{TRUE} it returns the negative 
#' spectral entropy. This is useful when maximizing forecastibility which is 
#' equivalent (up to an additive constant) to maximizing negative entropy. 
#' Default: \code{FALSE}.
#' @return 
#' \code{foreca.EM.h} returns non-negative real value (see References for details):
#' \itemize{
#'    \item entropy, if \code{weightvector.new = weightvector.current},
#'    \item an upper bound of that entropy for \code{weightvector.new},
#'          otherwise.
#' }
#' @seealso
#' \code{\link{weightvector2entropy_wcov}}
#' @export
#' @examples
#' 
#' foreca.EM.h(ww0, ff)       # iteration 0
#' foreca.EM.h(ww1, ff, ww0)  # min eigenvalue inequality
#' foreca.EM.h(ww1, ff)       # KL divergence inequality
#' one.step$value
#' 
#' # by definition of Omega, they should equal 1 (modulo rounding errors)
#' Omega(mvspectrum.output = f.ww0) / 100 + foreca.EM.h(ww0, ff)
#' Omega(mvspectrum.output = f.ww1) / 100 + foreca.EM.h(ww1, ff)
#' 

foreca.EM.h <- function(weightvector.new, f.U, 
                        weightvector.current = weightvector.new, 
                        f.current = NULL,
                        entropy.control = list(),
                        return.negative = FALSE) {
  
  # short as quadratic_form(apply(f.U*-log(weightvector.current), 2:3, sum),
  # weightvector.new)
  stopifnot(length(weightvector.current) == length(weightvector.new),
            all(f.current >= 0) || is.null(f.current))
  
  check_mvspectrum_normalized(f.U)
  if (is.null(f.current)) {
    f.current <- foreca.EM.E_step(f.U, weightvector.current)
  } else {
    check_mvspectrum_normalized(f.current)
  }

  num.series <- length(weightvector.new)  
  num.freqs <- length(f.current)

  entropy.control <- complete_entropy_control(entropy.control,
                                              num.outcomes = 2 * num.freqs)
  integrated.spectrum.entropy <- 
      weightvector2entropy_wcov(NULL, 
                                f.U = f.U,
                                f.current = f.current,
                                entropy.control = entropy.control)
  
  return(quadratic_form(integrated.spectrum.entropy, weightvector.new))
}

