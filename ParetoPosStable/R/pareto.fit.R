#' @title Fitting a Pareto distribution
#' @description It is an auxiliar function for fitting a Pareto distribution as a particular case of a Pareto Positive Stable distribution, allowing the scale parameter to be held fixed if desired.
#' @param x The vector of observations.
#' @param estim.method The estimation method, "MLE" or "OLS"
#' @param sigma The value of the scale parameter, if it is known; if the value is \code{NULL}, the parameter is estimated.
#' @param start Unused argument from \code{PPS.fit}.
#' @param \dots Other arguments.
#' 
#' @details This function is called by \code{PPS.fit()} when \code{Pareto} argument is \code{TRUE}.
#' 
#' @return A \code{PPSfit} Object.
#'
#' @references Sarabia, J.M and Prieto, F. (2009). The Pareto-positive stable distribution: A new descriptive model for city size data, \emph{Physica A: Statistical Mechanics and its Applications}, \bold{388}(19), 4179-4191.
#'
#' @seealso \code{\link{PPS.fit}}, \code{\link{coef.PPSfit}}, \code{\link{print.PPSfit}}, \code{\link{plot.PPSfit}}, \code{\link{GoF}}

#' @importFrom stats lm optimize

#' @export
pareto.fit <-
  function(x, estim.method = "MLE", sigma = NULL, start, ...){
    if (!is.element(estim.method,c("MLE", "OLS"))) stop("The estimation method in method.estim is not OLS nor MLE")
    
    if (missing(x) || length(x) == 0L || mode(x) != "numeric" || any(x<=0)) stop("'x' must be a non-empty numeric vector of positive values")
    if (any(!is.finite(x))) stop("'x' contains missing or infinite values")
    
    n <- length(x)
    Fn<-ecdf(x)
    y<-log(1-n*Fn(x)/(n+1))
    xName <- paste(deparse(substitute(x), 500), collapse = "\n")
    if (missing(start)) start <- NULL
    dots <- names(list(...))
    
    if (estim.method == "OLS"){
      if (is.null(sigma)){
        sce<-function(sigma){
          xx<-log(sigma)-log(x)
          ajuste <- lm(y ~ -1 + xx)
          sce<-sum(ajuste$residuals^2)
        }
        sigma <- optimize(f=sce,interval=c(0,min(x)))$minimum
        xx<-log(sigma)-log(x)
        ajuste <- lm(y ~ -1 + xx)
        lambda <- ajuste$coefficients[[1]]
        loglik <- n*(log(lambda)+lambda*log(sigma))-(lambda+1)*sum(log(x))
        return(structure(list(estimate = list(lambda = lambda, sigma = sigma), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS", Pareto = TRUE), class = "PPSfit"))
      }
      if (!is.null(sigma)){
        sc<-sigma
        xx<-log(sigma)-log(x)
        ajuste <- lm(y ~ -1 + xx)
        lambda <- ajuste$coefficients[[1]]
        loglik <- n*(log(lambda)+lambda*log(sigma))-(lambda+1)*sum(log(x))
        return(structure(list(estimate = list(lambda = lambda), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "OLS", sigma = sigma, Pareto = TRUE), class = "PPSfit"))
      }
    }
    
    if (estim.method == "MLE"){
      missing(start)
      if (is.null(sigma)){
        sigma <- min(x)
        lambda <- n/sum(log(x/sigma))
        loglik <- n*(log(lambda)+lambda*log(sigma))-(lambda+1)*sum(log(x))
        return(structure(list(estimate = list(lambda = lambda, sigma = sigma), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "MLE", Pareto = TRUE), class = "PPSfit"))
      }
      if (!is.null(sigma)){
        lambda <- n/sum(log(x/sigma))
        loglik <- n*(log(lambda)+lambda*log(sigma))-(lambda+1)*sum(log(x))
        return(structure(list(estimate = list(lambda = lambda), loglik = loglik, n = n, obs = x, obsName = xName, estim.method = "MLE", sigma = sigma, Pareto = TRUE), class = "PPSfit"))
      }
    }
  }