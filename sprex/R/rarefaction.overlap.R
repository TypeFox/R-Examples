#' @title Rarefaction Overlap
#' @description Calculate the percent of overlap between two 
#'   species estimate distributions where the larger sample size has been
#'   rarefied to match the smaller sample size.
#' 
#' @param x,y two vectors of species frequencies where the \code{i-th} element is the 
#'   number of species represented by only \code{i} samples.
#' @param f0.func function to use to calculate f0. Can be \code{\link{Chao1}},
#'   \code{\link{ACE}}, \code{\link{jack1}}, \code{\link{jack2}}, 
#'   \code{\link{iChao1}}, or \code{\link{Swor1}}.
#' @param n.rare sample size to rarefy both populations to. Must be <= the minimum
#'   sample size. If \code{NULL}, the minimum sample size is used. 
#' @param ... other arguments to \code{f0.func}.
#' 
#' @details Calculates the expected number of species and the standard 
#'   deviation for the smaller sample size of \code{x} and \code{y} using 
#'   the frequency distributions of each. The function then fits a gamma 
#'   distribution to each of these estimates, and returns the percent of overlap
#'   as the integral of the mininum value of the PDF for the two distributions. 
#'   Integration takes place from 0 to the largest quantile representing 
#'   0.99999 of either distribution.
#' 
#' @return a vector with the percent of overlap between the two distributions, 
#'   the sample size, and species estimates for the \code{x} and 
#'   \code{y} vectors.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @references Colwell, R.K., A. Chao, N.J. Gotelli, S.-Y. Lin, C.X. Mao, 
#'   R.L. Chazdon, and J.T. Longino. 2012. Models and estimators linking 
#'   individual-based and sample-based rarefaction, extrapolation and 
#'   comparison of assemblages. Journal of Plant Ecology 5(1):3-21.
#' 
#' @seealso \code{\link{discovery.curve}}
#' 
#' @examples
#' data(osa.old.growth)
#' data(osa.second.growth)
#' x <- expand.freqs(osa.old.growth)
#' y <- expand.freqs(osa.second.growth)
#' rarefaction.overlap(x, y, Chao1)
#' 
#' @importFrom stats dgamma qgamma integrate
#' @export
#' 
rarefaction.overlap <- function(x, y, f0.func, n.rare = NULL, ...) {
  # get minimum n
  n.min <- pmin.int(f.stats(x)["n"], f.stats(y)["n"])
  if(is.null(n.rare)) {
    n.rare <- n.min
  } else if(n.rare > n.min) {
    stop("'n.rare' must be <= minimum sample size in 'x' and 'y'")
  }
  
  # estimate number of species and standard deviation at minimum n
  s.ind.x <- expected.num.species(n.rare, x, f0.func, ...)
  s.ind.y <- expected.num.species(n.rare, y, f0.func, ...)
  mu.xy <- unname(c(s.ind.x["s.ind"], s.ind.y["s.ind"]))
  sd.xy <- unname(c(s.ind.x["sd.s.ind"], s.ind.y["sd.s.ind"]))
  
  # return NA if can't estimate
  result <- c(
    n.rare = n.rare, x.est = mu.xy[1], x.sd = sd.xy[1], 
    y.est = mu.xy[2], y.sd = sd.xy[2]
  )
  if(any(is.na(result)) | any(sd.xy == 0)) return(c(pct.overlap = NA, result))

  # convert mu and sd to gamma parameters
  sh <- (mu.xy / sd.xy) ^ 2
  sc <- (sd.xy ^ 2) / mu.xy
  # integrate minimum PDF
  min.pdf <- function(x, shape, scale) {
    pmin(dgamma(x, shape = shape[1], scale = scale[1]), 
         dgamma(x, shape = shape[2], scale = scale[2])
    )
  }
  upper <- max(qgamma(0.99999, shape = sh, scale = sc))
  
  c(pct.overlap = integrate(min.pdf, 0, upper, shape = sh, scale = sc)$value,
    result
  )
}