#' @title Bootstrap Assemblage of Species
#' @description Create bootstrap assemblage of species.
#' 
#' @param f a vector of species frequencies where \code{f[i]} is the number of 
#'   species represented by only \code{i} samples.
#' @param f0.func function calculating the unobserved number of 
#'   species (\code{f0}).
#' @param n.boot number of bootstrap replicates.
#' @param ... other arguments to \code{f0.func}.
#' 
#' @return a list of bootstrap replicates of species frequencies.
#' 
#' @author Eric Archer \email{eric.archer@@noaa.gov} 
#' 
#' @references Chao, A., N.J. Gotelli, T.C. Hsieh, E.L. Sander, K.H. Ma, 
#'   R.K. Colwell, and A.M. Ellison. 2014. Rarefaction and extrapolation with
#'   Hill numbers: a framework for sampling and estimation in species 
#'   diversity studies. Ecological Monographs 84(1):45-67.
#'   
#' @importFrom stats rmultinom   
#' @export

bootstrap.assemblage <- function(f, f0.func, n.boot = 500, ...) {
  if(length(f) == 1) f <- c(f, 0)
  f0.est <- f0.func(f, ...)
  n <- unname(f0.est["n"])
  f0 <- unname(f0.est["f0"])
  f0.star <- ceiling(f0)

  # calcuate Cind(n)
  term.1 <- f[1] / n
  term.2 <- (n - 1) * ifelse(f[2] > 0, f[1], f[1] - 1)
  term.3 <- term.2 + 2 * ifelse(f[2] > 0, f[2], 1)
  Cind.n <- 1 - term.1 * term.2 / term.3

  # calculate tuned relative abundance for each species
  sample.freq <- species.to.sample.freq(f)
  sample.freq <- sample.freq[sample.freq > 0]
  term.1 <- sample.freq / n
  term.2 <- (1 - term.1) ^ n
  term.3 <- sum(term.1 * term.2)
  lambda.hat <- (1 - Cind.n) / term.3
  est.prob.i <- term.1 * (1 - (lambda.hat * term.2))
  est.prob.unseen <- rep((1 - Cind.n) / f0.star, f0.star)

  # generate bootstrap assemblage
  probs <- c(est.prob.i, est.prob.unseen)
  if(any(is.na(probs) | is.nan(probs))) return(NULL)
  boot.sample <- rmultinom(n.boot, n, probs)
  # convert bootstrap assemblage of species frequencies to list of f_i (number of species represented by i individuals)
  lapply(1:ncol(boot.sample), function(i) sample.to.species.freq(boot.sample[, i]))
}
