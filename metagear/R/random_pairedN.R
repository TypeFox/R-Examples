#' Random generation of paired sample sizes (N) for study outcomes.  
#'
#' Generates random paired sample sizes (N).  For example, sample sizes for a
#'    treatment group and samples sizes for a control group.  These paired N 
#'    are often correlated within studies.
#'
#' @param K Number of paired sample sizes to generate.
#' @param mean The lambda (dispersion parameter) of a Poisson distribution.  The
#'    default is 15, which will generate sample sizes that on average will center
#'    around N = 15.
#' @param min A non-negative integer that specifies the minimum sample size that
#'    can be generated.  Default is N = 3. 
#' @param correlation A correlation ranging from zero to one that specifies how
#'    'similar' the paired sample sizes will be to one another.  Default is 0.95 
#'    (i.e. the paired sample sizes will be highly correlated).
#'
#' @return A data table of paired random sample sizes (N).
#'
#' @importFrom stats pnorm qpois
#' @export random_pairedN

random_pairedN <- function(K, 
                           mean = 15,
                           min = 3,
                           correlation = 0.95) {
  
  corMatrix <- matrix(c(1.0, correlation, correlation, 1.0), ncol = 2, nrow = 2)
  lambda <- matrix(c(mean, mean), ncol = 1, nrow = 2)
  thek <- K
  theN <- c()
  
  while(thek != 0) {
    theN <- rbind(GenerateMultivariatePoisson(2, thek, corMatrix, lambda), theN)
    theN <- theN[theN[, 1] >= min & theN[, 2] >= min, ]
    thek <- K - nrow(theN)
  }
  
  return(theN)
}
