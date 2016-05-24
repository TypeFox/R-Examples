#' Random Generation for MoTBFs
#' 
#' Random generation for mixtures of truncated basis functions defined in a specific domain.
#' The inversion method is used. It is a technique to get random samples from a 
#' probability distribution.
#' 
#' @name MoTBF-Distribution
#' @rdname MoTBF-Distribution
#' @param size A \code{"numeric"} value indicating the number of records to generate.
#' @param fx An object of class \code{"motbf"}.
#' @param domain A \code{"numeric"} vector with two values indicating the range where generating
#' the data sample. By default it is \code{NULL} and the range is taken of the function, \code{fx}.
#' @param data A \code{"numeric"} vector which contains the data which we want to compare with the generated 
#' function. By default it is \code{NULL}, if not the empirical cumulative distribution is plotted and the 
#' Kolmogorov Smirnov test is used to compare the generated sample and the data.
#' @return \code{rMoTBF()} returns a \code{"numeric"} vector containing the simulated values. \code{inversionMethod()} 
#' returns a list with the simulated values and the results of the test, it also shows a plot with the \bold{cdf}
#' of the original data and the generated one by screen.
#' @seealso \link{integralMoTBF}
#' @examples
#' 
#' ## 1. EXAMPLE
#' ## Data
#' X <- rnorm(1000, mean = 5, sd = 3)
#' 
#' ## Learning
#' f <- univMoTBF(X, POTENTIAL_TYPE="MOP", nparam=10)
#' plot(f, xlim = f$Domain)
#' 
#' ## Random sample
#' Y <- rMoTBF(size = 500, fx = f)
#' ks.test(X,Y)
#' 
#' ## Plots
#' hist(Y, prob = TRUE, add = TRUE)
#' 
#' ## 2. EXAMPLE 
#' ## Data
#' X <- rweibull(5000, shape=2)
#' 
#' ## Learning
#' f <- univMoTBF(X, POTENTIAL_TYPE="MOP", nparam=10)
#' plot(f, xlim = f$Domain)
#' 
#' ## Random sample
#' inv <- inversionMethod(size = 500, fx = f, data = X)
#' attributes(inv)
#' inv$test
#' Y <- inv$sample 
#' 
#' ## Plots
#' plot(f, xlim = f$Domain)
#' hist(Y, prob = TRUE, add = TRUE)
#' 

#' @rdname MoTBF-Distribution
#' @export
rMoTBF <- function(size, fx, domain = NULL)
{
  if(is.null(domain)) domain <- fx$Domain
  if(is.null(domain)) return(cat("Domain is required."))
  CDF <- integralMoTBF(fx)
  intmin <- as.function(CDF)(min(domain))
  mUnif <- runif(size)
  simulatedValues <- sampleMoTBF <- sapply(1:size, function(i)
                     uniroot(as.function(motbf(paste(CDF,ifelse
                     ((-intmin-mUnif[i])>=0, "+", ""),-intmin-
                     mUnif[i], sep=""))), range(domain))$root)
  return(simulatedValues)
}

#' @name MoTBF-Distribution
#' @export
inversionMethod <- function(size, fx, domain = NULL, data = NULL)
{
  simulatedValues <- rMoTBF(size, fx, domain)
  
  plot(ecdf(simulatedValues), cex = 0, main = "")
  if(!is.null(data)){
    plot(ecdf(data), col="red", cex = 0, main = "", add = T)
    par(xpd = TRUE)
    legend(0, 1.4, c(expression(F(X)),expression(F(Simulated_Values))), 
         cex=0.8, col = c("red", "black"), lty = c(1,1), lwd = c(1,1), 
         bty = "n")
    test <- ks.test(data, simulatedValues)
    simulatedValues <- list(sample = simulatedValues, test = test)
  }
  
  return(simulatedValues)
}

