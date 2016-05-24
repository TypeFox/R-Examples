#' @title Bootstrapping for cdf quantile regression
#' @aliases qrBoot
#' @description \code{qrBoot} provides a simple bootstrapping method for estimating the parameters of a cdf quantile regression model.
#' 
#' @param object The fitted cdfqr model object
#' @param rn The sample size of bootstrap samples
#' @param f A function whose one argument is the name of a cdfqr object that will be applied to the updated cdfqr object to compute the statistics of interest. The default is coef.
#' @param R Number of bootstrap samples. 
#' @param ci The confidence interval level to obtain the boostrap confidence intervals
#'
#' @return A matrix that includes the original statistics, boostrap means, and boostrap confidence intervals
#' @export
#'
#' @examples
#' data(cdfqrExampleData)
#'fit <- cdfquantreg(crc99 ~ vert | confl, 't2', 't2', data = JurorData)
#'qrBoot(fit, rn = 50, R = 50)
#' 
qrBoot <- function (object, rn,  f = coef, R = 500, ci = 0.95) 
{
  # Extract the call from the original model
  call <- object$call
  dat <- eval(call$data)
  fd1 <- object$family$fd
  sd1 <- object$family$sd
  n <- nrow(dat)
  stats0 <- do.call(f, list(object)) # The original stats
  

  # Determin if seed has been set
  if (!exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)) 
    runif(1)
  seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
  
 
  # Initiate the bootstrap results 
  stats <- NULL
  
  for (i in 1:R) {
    #get bootstrap sample
    ind <- sample(1:n, size = rn, replace = TRUE)
    dat1 <- as.data.frame(dat[ind, ])
    
    #modify the call for cdfqr function with the new dataset
    mod <- update(object, .~., fd=fd1, sd=sd1, data = dat1)

    #extract results
    stats_temp <- do.call(f, list(mod))  
    
    stats <- rbind(stats, stats_temp)
  }
  
 
  boot.mean <- colMeans(stats)
  boot.ci <- apply(stats, 2, quantile, probs = c((1-ci)/2, 1-(1-ci)/2)) 
  
  boot.output <- data.frame(original = stats0, 
                            boot.mean = boot.mean,
                            ci.low = boot.ci[1,],
                            ci.upper = boot.ci[2,])
  
  boot.output
}