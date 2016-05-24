#' Accurate Computation
#' @description Computes the accurate criterion of smoothed (fitted) values.
#' @param x a numeric vector of original values.
#' @param x.hat a numeric vector of smoothed (fitted) values.
#' @param k the number of parameters in obtaining the smoothed (fitted) values.
#' @param output a logical value indicating to print the results in R console. The default is
#' \code{TRUE}.
#' @details See \url{http://www.dms.umontreal.ca/~duchesne/chap12.pdf} in page 616 - 617 for
#' the details of calculations for each criterion.
#' @note If the model fits the series badly, the model error sum of squares \code{SSE}
#' may be larger than \code{SST} and the \code{R.squared} or \code{RW.R.squared} statistics 
#' will be negative. The \code{RW.R.squared} uses the random walk model for the purpose of 
#' comparison.
#' @return A vector containing the following components:
#' \item{SST}{the total sum of squares.}
#' \item{SSE}{the sum of the squared residuals.}
#' \item{MSE}{the mean squared error.}
#' \item{RMSE}{the root mean square error.}
#' \item{MAPE}{the mean absolute percent error.}
#' \item{MPE}{the mean percent error.}
#' \item{MAE}{the mean absolute error.}
#' \item{ME}{the mean error.}
#' \item{R.squared}{R^2 = 1 - SSE/SST.}
#' \item{R.adj.squared}{the adjusted R^2.}
#' \item{RW.R.squared}{the random walk R^2.}
#' \item{AIC}{the Akaike's information criterion.}
#' \item{SBC}{the Schwarz's Bayesian criterion.}
#' \item{APC}{the Amemiya's prediction criterion}
#' @author Debin Qiu
#' @examples X <- matrix(rnorm(200),100,2)
#' y <- 0.1*X[,1] + 2*X[,2] + rnorm(100)
#' y.hat <- fitted(lm(y ~ X))
#' accurate(y,y.hat,2)
#' @importFrom stats embed
#' @export
accurate <- function(x,x.hat,k,output = TRUE)
{
  n <- length(x)
  SST <- sum((x - mean(x))^2)
  SSE <- sum((x - x.hat)^2)
  MSE <- SSE/(n - k)
  RMSE <- sqrt(MSE)
  MAPE <- 100*mean(abs((x - x.hat)/x))
  MPE <- 100*mean((x - x.hat)/x)
  MAE <- mean(abs(x - x.hat))
  ME <- mean(x - x.hat)
  R2 <- 1 - SSE/SST
  ADJ.R2 <- 1 - (n - 1)*(1 - R2)/(n - k)
  z <- embed(x,2)
  RW.R2 <- 1 - (n - 1)*SSE/(n*sum(z[,1] - z[,2] - mean(z[,1] - z[,2])))
  AIC <- n*log(SSE/n) + 2*k
  SBC <- n*log(SSE/n) + k*log(n)
  APC <- ((n + k)/(n*(n - k)))*SSE
  result <- c(SST,SSE,MSE,RMSE,MAPE,MPE,MAE,ME,R2,ADJ.R2,RW.R2,AIC,SBC,APC)
  names(result) <- c("SST","SSE","MSE","RMSE","MAPE","MPE","MAE","ME","R.squared",
                     "R.adj.squared","RW.R.squared","AIC","SBC","APC")
  if (output) {
    cat("SST:",SST, "; SSE:", SSE, "; MSE:", MSE, "; RMSE:", RMSE, "\n")
    cat("MAPE:", MAPE, "; MPE:", MPE, "; MAE:", MAE, "; ME:", ME, "\n")
    cat("R.squared:", R2, "; R.adj.squared:", ADJ.R2, "; RW.R.squared:", RW.R2,"\n")
    cat("AIC:", AIC, "; SBC:", SBC, "; APC:", APC)
  }
  accurate <- result
}