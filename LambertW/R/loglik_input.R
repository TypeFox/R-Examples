#' @rdname loglik-LambertW-utils
#' @description
#' \code{loglik_input} computes the log-likelihood of various distributions for 
#' the parameter \eqn{\boldsymbol \beta} given the data \code{x}. This can be 
#' used independently of the Lambert W x F framework to compute 
#' the log-likelihood of parameters for common distributions.
#' 
#' @inheritParams common-arguments
#' @param dX optional; density function of \code{x}. Common distributions are
#'     already built-in (see \code{distname}). If you want to supply your own
#'     density, you \strong{must} supply a function of \code{(x, beta)} and set
#'     \code{distname = "user"}.
#' @param log.dX optional; a function that returns the logarithm of the density
#'     function of \code{x}. Often -- in particular for exponential families --
#'     the \eqn{\log} of \eqn{f_X(x)} has a simpler form (and is thus faster to
#'     evaluate).
#' @param x a numeric vector of real values (the \emph{input} data).
#' @export
loglik_input <- function(beta, x, distname, dX = NULL, 
                         log.dX = function(x, beta) log(dX(x, beta))) {
  
  if (is.null(dX) && is.null(log.dX))  {
    stop("Please specify either the density function 'dX = ...' or \n",
         " (preferably) its logarithnm '.log_dX = ...'. \n ", 
         " In the form: 'dX = function(x) log(mydensity(x, params = beta))', ",
          "where beta is the parameter vector of 'mydensity' and specified as another ", 
          "argument of 'loglik_input'.")
  }
  if (distname != "user") {
    check_distname(distname)
    check_beta(beta, distname)
    names(beta) <- get_beta_names(distname)
  } else {
    .log_dX <- log.dX
  }
  switch(distname,
         cauchy = {
           .log_dX <- function(xx, beta = beta) {
             dcauchy(xx, location = beta[1], scale = beta[2], log = TRUE)
           }
         },
         chisq = {
           .log_dX <- function(xx, beta = beta) {
             dchisq(xx, df = beta[1], log = TRUE)
           }
         },
         exp = {
           .log_dX <- function(xx, beta = beta) {
             dexp(xx, rate = beta[1], log = TRUE)
           }
         },
         gamma = {
           .log_dX <- function(xx, beta = beta) {
             dgamma(xx, shape = beta["shape"], scale = beta["scale"], log = TRUE)
           }
         },
         normal = {
           .log_dX <- function(xx, beta = beta) {
             dnorm(xx, mean = beta[1], sd = beta[2], log = TRUE)
           }
         },
         t = {
           .log_dX <- function(xx, beta = beta) {
             dt((xx - beta["location"])/beta["scale"], df = beta["df"], log = TRUE) - log(beta["scale"])
           }
         },
         unif = {
           .log_dX <- function(xx, beta = beta) {
             dunif(xx, min = beta[1], max = beta[2], log = TRUE)
           }
         },
         user = {
     
         }
         )
   
  loglik <- sum(.log_dX(x, beta = beta))
  return(loglik)
} 
