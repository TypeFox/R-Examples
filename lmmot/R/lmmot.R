##' Fit right censored Multiple Ordinal Tobit (MOT) Model.
##'
##' Fit right censored Multiple Ordinal Tobit (MOT) model. 
##' The model is a right censored Tobit model with multiple ordinal categories for latent values above the threshold, the threshold is therefore replaced by a threshold vector.
##' 
##' For the latent variable a linear model with independent and identically distributed non-systematic and homoscedastic errors is assumed.
##' 
##' If the threshold is of length 1, the model is equivalent to the standard right censored Tobit model. 
##' 
##' The data is fitted with the Maximum Likelihood method.
##' 
##' @title Multiple Ordinal Tobit Model
##' 
##' @param formula Object of class formula describing the model.
##' @param data Optional data frame or environment containing the variables in the model.
##' @param threshold Vector of thresholds in the model.
##' @param stdEr Method for standard error estimation. 
##' Use "fisher" for estimation using the inverse of the Fisher information matrix or "hessian" for estimation using the Hessian matrix. 
##' @param ... Further arguments passed to the maximum likelihood estimation function \link[maxLik]{maxLik}.  
##'
##' @return lmmot object: \link[maxLik]{maxLik} object with additional fields: 
##' \itemize{
##'  \item{censoring: Number of obeservations in the censoring intervals.}
##'  \item{fisherInfo: Fisher information matrix.}
##'  \item{stdEr: Standard errors for estimated coefficients.}
##'  \item{tval: Value for t statistic in Wald test.}
##'  \item{pval: p-value in Wald test.}
##'  \item{fitted.values: Fitted values of the estimated model.}
##'  \item{residuals: Residuals of the estimated model.}
##' }
##' 
##' @export
##' @import stats maxLik MASS
##' @seealso \link[stats]{lm} \link[maxLik]{maxLik}
##' @examples
##' # Random data for x
##' N <- 100
##' x <- rnorm(N, 25, 10)
##' 
##' # Simulate data for latent variable ystar with simple linear model
##' beta_0 <- 60
##' beta_1 <- 1
##' sigma <- 8
##' ystar <- beta_0 + beta_1*x + rnorm(N, 0, sigma)
##' 
##' # Simulate censoring for observed variable y
##' y <- ystar
##' y[y >= 100] <- 100
##' y[(y >= 90) & (y < 100)] <- 90
##' y[(y >= 80) & (y < 90)] <- 80
##' 
##' # MOT regression with observed variable y
##' mot.fit <- lmmot(y ~ x, threshold = c(80, 90, 100))
##' 
##' # Show details
##' summary(mot.fit)
##' 
##' # Compare real data with model fit 
##' plot(x, ystar)
##' abline(coefficients(mot.fit)[1:2])
##' @author Marvin N. Wright

lmmot <- function(formula, data=sys.frame(sys.parent()), threshold, stdEr = "fisher", ...) {
  
  # check arguments
  if (class(formula) != "formula") {
    stop("invalid formula.")
  }
  if (!is.numeric(threshold)) {
    stop("threshold must be numeric.")
  }
  
  # get x and y
  y <- model.frame(formula,data)[,1]
  x <- model.matrix(formula,data)
  n <- length(y) 
  p <- ncol(x)
  
  # count observations in intervals
  K <- length(threshold)
  tau <- c(-Inf,threshold,Inf)
  n.cens <- vector("numeric",K)
  for (k in 1:(K+1)) {
    n.cens[k] <- sum((y >= tau[k])*(y < tau[k+1]))
  }
  cens.names <- paste("c", 1:K, sep="")
  
  # use normal linear model for starting values
  fit <- lm(formula,data)
  beta <- coef(fit)
  sigma <- sqrt(1/n * sum(residuals(fit)^2))
  start <- c(beta,sigma)
  
  # compute ml estimator
  result <- maxLik(logLik=motLogLik, grad=motGradient, hess=motHessian,
                   xx=t(x), y=y, tau=threshold, start=start, ...) 
  
  # compute fisher info
  fish <- motFisher(result$estimate, t(x), threshold)
  
  # compute standard error 
  if (stdEr == "fisher") {
    stdEr <- sqrt(diag(ginv(fish)))
  } else if (stdEr == "hessian") {
    stdEr <- stdEr(result)
  } else {
    stop("unknown method for stdEr.")
  }
          
  # perform wald-test
  tval <- result$estimate / stdEr
  pval <- 2 * pt(-abs(tval), df=n-p)
    
  # compute fitted values and residuals
  fitted <- x %*% result$estimate[-length(result$estimate)]
  residuals <- y - fitted
    
  # set result elements
  class(result) <- c("lmmot", class(result))
  result$call <- match.call()
  names(result$estimate)[length(names(result$estimate))] <- "sigma"
  result$censoring <- c(n.cens, n)
  names(result$censoring) <- c("uncensored", cens.names, "total")
  result$fisherInfo <- fish
  result$stdEr <- stdEr
  result$tval <- tval
  result$pval <- pval
  result$fitted.values <- fitted
  result$residuals <- residuals
  
  return(result)
}
