#'MLE for Gamma distribution
#'
#'Maximum likelihood estimation of Gamma distributed observations
#'distribution parameters
#'
#'@param g a list of Gamma distributed observation.
#'
#'@importFrom stats uniroot
#'
#'@export
#'
#'@examples
#'
#' g_list <- list()
#' for(i in 1:1000){
#'  g_list <- c(g_list, rgamma(1, shape=100, rate=5))
#' }
#'
#' mle <- MLE_gamma(g_list)
#' mle
#'
MLE_gamma <- function(g){
  N <- length(g)

  a_mle <- try(stats::uniroot(function(a){(N*mean(log(g))
                                           - N*digamma(a)
                                           - N*log(mean(g))
                                           + N*log(a)
  )}, lower = 0.000001, upper=1E9)$root, TRUE)
  if(inherits(a_mle, "try-error")){
    #either
    a_mle <- 0.0001
    b_mle <- 0.0001
    warning("Unable to estimate Gamma hyperpriors properly
            (this can happen when only a few clusters are fitted).
            => Non informative values are returned instead")
  }else{
    b_mle <- mean(g)/a_mle
  }

  return(list("shape"=a_mle, "scale"=b_mle, "rate"=1/b_mle))
}
