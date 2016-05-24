#' Random generation of sample sizes (N) for study outcomes.  
#'
#' Generates random sample sizes (N) by either sampling from a Negative Binomial 
#'    or Poisson distribution.       
#'
#' @param K Number of sample sizes to generate.
#' @param method A string that defines what sampling distribution to generate
#'    random N.  The default is \code{"NegativeBinomial"} but a \code{"Poisson"} 
#'    distribution can also be used.
#' @param mean The population mean (mu) if \code{"NegativeBinomial"}, or the lambda 
#'    (dispersion parameter) if \code{"Poisson"}.  The default is 15, which will 
#'    generate sample sizes that on average will center around N = 15.
#' @param min A non-negative integer that specifies the minimum sample size that
#'    can be generated.  Default is N = 3. 
#' @param NB_size Dispersion parameter for the \code{"Negative Binomial"} distribution
#'    that must be strictly positive, but need not be integer.  Default is 15, 
#'    which creates a long tail for random N's ranging to about N = 60.  Increase
#'    value to create a longer tail of random sample sizes.
#'
#' @return A vector of random sample sizes (N).
#'
#' @importFrom stats rnbinom rpois
#' @export random_N

random_N <- function(K, 
                     method = "NegativeBinomial",
                     mean = 15,
                     min = 3,
                     NB_size = 15) {
  
  thek <- K
  theN <- c()
  
  while(length(theN) != K) {
    if(method == "NegativeBinomial") 
      theN <- c(rnbinom(thek, mu = mean, size = NB_size), theN)
    else if (method == "Poisson") theN <- c(rpois(thek, lambda = mean), theN)
    else message("error") #need to update
    theN <- theN[theN >= min]
    thek <- K - length(theN)
  }
  
  return(theN) 
}