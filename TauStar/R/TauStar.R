#' Efficient Computation and Testing of the t* Statistic of Bergsma and Dassios
#'
#' Computes the t* statistic corresponding to the tau star population
#' coefficient introduced by Bergsma and Dassios (2014) <DOI:10.3150/13-BEJ514>
#' and does so in  O(n^2*log(n)) time following the algorithm of Weihs,
#' Drton, and Leung (2016) <DOI:10.1007/s00180-015-0639-x>. Also allows for
#' independence testing using the asymptotic distribution of t* as described by
#' Nandy, Weihs, and Drton (2016) <http://arxiv.org/abs/1602.04387>.
#' To directly compute the t* statistic see the function tStar. If otherwise
#' interested in performing tests of independence then see the function
#' tauStarTest.
#'
#' @references
#' Bergsma, Wicher; Dassios, Angelos. A consistent test of independence based
#' on a sign covariance related to Kendall's tau. \emph{Bernoulli} 20 (2014), no.
#' 2, 1006--1028.
#' \cr\cr
#' Luca Weihs, Mathias Drton, and Dennis Leung. Efficient Computation of the
#' Bergsma-Dassios Sign Covariance. \emph{Computational Statistics}, x:x-x,
#' 2016. to appear.
#' \cr\cr
#' Preetam Nandy, Luca Weihs, and Mathias Drton. Large-Sample Theory for the
#' Bergsma-Dassios Sign Covariance. arXiv preprint arXiv:1602.04387. 2016.
#'
#' @examples
#' \dontrun{
#' library(TauStar)
#'
#' # Compute t* for a concordant quadruple
#' tStar(c(1,2,3,4), c(1,2,3,4)) # == 2/3
#'
#' # Compute t* for a discordant quadruple
#' tStar(c(1,2,3,4), c(1,-1,1,-1)) # == -1/3
#'
#' # Compute t* on random normal iid normal data
#' set.seed(23421)
#' tStar(rnorm(4000), rnorm(4000)) # near 0
#'
#' # Compute t* as a v-statistic
#' set.seed(923)
#' tStar(rnorm(100), rnorm(100), vStatistic=TRUE)
#'
#' # Compute an approximation of tau* via resampling
#' set.seed(9492)
#' tStar(rnorm(10000), rnorm(10000),
#'       resample=TRUE, sampleSize=30, numResamples=5000)
#'
#' # Perform a test of independence using continuous data
#' set.seed(123)
#' x = rnorm(100)
#' y = rnorm(100)
#' testResults = tauStarTest(x,y)
#' print(testResults$pVal) # big p-value
#'
#' # Now make x and y correlated so we expect a small p-value
#' y = y + x
#' testResults = tauStarTest(x,y)
#' print(testResults$pVal) # small p-value
#' }
#'
#' @importFrom Rcpp evalCpp
#' @useDynLib TauStar
#' @import stats
"_PACKAGE"
#> [1] "_PACKAGE"
