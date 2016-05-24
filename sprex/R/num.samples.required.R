#' @title Number of Samples Required
#' @description Calculate the additional number of samples to required to observe a 
#'   given proportion of the total number of species.
#' 
#' @param g propotion of total number of species.
#' @param f a vector of species frequencies where \code{f[i]} is the number of 
#'   species represented by only \code{i} samples.
#' @param f0.func a function that computes the number of unobserved 
#'   species (f0).
#' @param ... other arguments to \code{f0.func}.
#' 
#' @return a vector containing of the estimated additional number of samples 
#'   (\code{m.g}) required to observe \code{g} percent of the total number 
#'   of species.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @references Eqn 12 in Chao, A., R.K. Colwell, C.-W. Lin, and N.J. Gotelli. 2009. 
#'   Sufficient sampling for asymptotic minimum species richness estimators. 
#'   Ecology 90(4):1125-1133. \cr
#'   Eqn 11 in Colwell, R.K., A. Chao, N.J. Gotelli, S.-Y. Lin, 
#'   C.X. Mao, R.L. Chazdon, and J.T. Longino. 2012. Models and estimators 
#'   linking individual-based and sample-based rarefaction, extrapolation and 
#'   comparison of assemblages. Journal of Plant Ecology 5(1):3-21.
#'    
#'  
#' @examples
#' data(osa.old.growth)
#' f <- expand.freqs(osa.old.growth)
#' num.samples.required(0.6, f = f, f0.func = Chao1)
#' 
#' @importFrom stats optim
#' @importFrom swfscMisc isBetween
#' @export

num.samples.required <- function(g, f, f0.func, ...) {
  if(!isBetween(g, 0, 1)) stop("'g' must be between 0 and 1")
  x <- f0.func(f, ...)
  s.est <- unname(x["s.est"])
  
  # calculate additional m individuals required to detect g proportion of 
  #   s.est (Eqn. 11, Colwell et al 2012)
  m.g <- if(x["f0"] > 0) {
    if(f[2] > 0) {
      term.1 <- (x["n"] * f[1]) / (2 * f[2])
      term.2 <- x["f0"] / ((1 - g) * s.est)
      term.1 * log(term.2)
    } else {
      warning("since f2 == 0, number of samples based on optimizaton of Colwell et al 2012 Eqn. 9")
      result <- optim(
        par = c(m.star = 1),
        fn = function(m.star, s.g, f0, f1, n, s.obs) {
          s.ind <- .s.ind.n.m(f0, f1, n, m.star, s.obs)
          unname(abs(s.ind - s.g))
        },
        method = "L-BFGS-B", lower = 0,
        s.g = g * s.est, f0 = x["f0"], f1 = f[1], n = x["n"], s.obs = x["s.obs"]
      )
      if(result$convergence != 0) {
        msg <- paste("equation 9 optimization did not converge, with code ",
                     result$convergence, ": ", result$message)
        warning(msg)
        NA
      } else result$par
    }
  } else 0
  c(m.g = unname(m.g), g = g, x)
}
