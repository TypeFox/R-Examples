#' @title MedCouple Estimator
#' 
#' @description
#' A robust measure of asymmetry. See References for details.
#' 
#' @param x numeric vector; if length > 3,000, it uses a random subsample
#'     (otherwise it takes too long to compute as calculations are of order
#'     \eqn{N^2}.)
#' @param seed numeric; seed used for sampling (when \code{length(x) >
#'     3000}).
#' @return float; measures the degree of asymmetry
#' @seealso \code{\link{test_symmetry}}
#' @references Brys, G., M. Hubert, and A. Struyf (2004). \dQuote{A robust
#'     measure of skewness}. Journal of Computational and Graphical Statistics
#'     13 (4), 996 - 1017.
#' @keywords univar
#' @export
#' @examples
#' 
#' # a simulation
#' kNumSim <- 100
#' kNumObs <- 200
#' 
#' ################# Gaussian (Symmetric) #### 
#' A <- t(replicate(kNumSim, {xx <- rnorm(kNumObs); c(skewness(xx), medcouple_estimator(xx))}))
#' ########### skewed LambertW x Gaussian #### 
#' tau.s <- gamma_01(0.2) # zero mean, unit variance, but positive skewness
#' rbind(mLambertW(theta = list(beta = tau.s[c("mu_x", "sigma_x")], 
#'                              gamma = tau.s["gamma"]), 
#'                 distname="normal"))
#' B <- t(replicate(kNumSim, 
#'                  {
#'                    xx <- rLambertW(n = kNumObs, 
#'                                    theta = list(beta = tau.s[c("mu_x", "sigma_x")], 
#'                                                 gamma = tau.s["gamma"]), 
#'                                    distname="normal")
#'                    c(skewness(xx), medcouple_estimator(xx))
#'                  }))
#'                   
#' colnames(A) <- colnames(B) <- c("MedCouple", "Pearson Skewness")
#' 
#' layout(matrix(1:4, ncol = 2))
#' plot(A, main = "Gaussian")
#' boxplot(A)
#' abline(h = 0)
#' 
#' plot(B, main = "Skewed Lambert W x Gaussian")
#' boxplot(B)
#' abline(h = mLambertW(theta = list(beta = tau.s[c("mu_x", "sigma_x")], 
#'                                   gamma = tau.s["gamma"]), 
#'                      distname="normal")["skewness"])
#' 
#' colMeans(A)
#' apply(A, 2, sd)
#' 
#' colMeans(B)
#' apply(B, 2, sd)
#' 

medcouple_estimator <- function(x, seed = sample.int(1e6, 1)) {

  stopifnot(is.numeric(x),
            is.numeric(seed),
            length(seed) == 1,
            seed > 0)
  set.seed(seed)
  if (length(x) > 3000) {
    warning("medcouple_estimator() is too slow for samples larger than 3000. ",
            "Using a subsample of 3000 instead.")
    x <- sample(x, size = 3000)
  }
  #### kernel
  .kernel <- function(z1, z2) {
    (z2 + z1)/(z2 - z1)
  }
  
  z <- sort(x, decreasing = TRUE)
  z <- z - median(z)
  
  Z.neg <- z[z < 0]
  Z.pos <- z[z > 0]

  return(median(outer(Z.neg, Z.pos, .kernel)))
} 
