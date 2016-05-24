#' Random generation of correlation coefficients.  
#'
#' Generates random correlation coefficients (r or Pearson product-moment 
#' correlation coefficients) and their variances (Pearson 1895).  Also provides 
#' Fisher z-transformed correlation coefficients (Fisher 1915).       
#'
#' @param K Number of effect sizes to generate.
#' @param correlation The mean population correlation coefficient (rho) to
#'    simulate.  Must range between -1 to 1.
#' @param N The number of samples used to estimate each correlation coefficient.
#'    When a non-negative integer, all r will be estimated using the same N.  A
#'    vector of unequal N's can also be taken; if so, K will be ignored and the 
#'    number of randomly generated r will equal the length of that vector.  
#' @param Fisher_Z When \code{TRUE}, also returns the Fisher z-transformed 
#'    correlation coefficients and their variances (Fisher 1915).
#'
#' @return A data table with columns of random effect sizes (r), their variances
#'    and sample sizes.
#'
#' @examples
#'    random_r(K = 5, correlation = 0.5, N = 50)
#'
#' @references Pearson, K. 1895. Notes on regression and inheritance in the 
#'    case of two parents. Proceedings of the Royal Society of 
#'    London 58: 240-242.
#' @references Fisher, R.A. 1915. Frequency distribution of the values of the
#'    correlation coefficient in samples of an indefinitely large 
#'    population. Biometrika 10: 507-521.
#'
#' @importFrom MASS mvrnorm
#' @importFrom stats cor
#' @export random_r

random_r <- function(K = 100,
                     correlation = 0.5,
                     N = 10, 
                     Fisher_Z = FALSE) {
  
  #add warning message N < 3
  
  if(length(N) != 1) {
    samples <- sum(N)
    theGroups <- rep(1:length(N), N)
    finalN <- N
    #add warning message that K is ignored
  }
  else {
    samples <- K * N
    theGroups <- rep(1:K, each = N)
    finalN <- N
  }
  
  C <- matrix(c(1, correlation, correlation, 1), nrow = 2, ncol = 2)  
  random_XY <- data.frame(mvrnorm(n = samples, mu = c(0, 0), C), 
                          group = theGroups)  
  theCor <- unlist(lapply(split(random_XY, random_XY$group), 
                                      function(x) return(cor(x[,1], x[,2]))))
  random_r <- data.frame(r = theCor, 
                         var_r = ((1.0 - theCor ^ 2.0) ^ 2.0)/(finalN - 2.0), 
                         N = finalN)
  
  if(Fisher_Z == TRUE) {
    random_r <- cbind(random_r, 
                      data.frame(Z = 0.5 * log((1.0 + random_r[, 1])/(1.0 - random_r[, 1])),
                                 var_Z = 1.0 / (random_r[, 3] - 3)))
  }
  
  return(random_r)
}
