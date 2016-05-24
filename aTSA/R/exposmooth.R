#' Simple Exponential Smoothing
#' @description Performs a simple exponential smoothing for univariate time series
#' with no trend or seasonal pattern.
#' @param x a numeric vector or univariate time series.
#' @param trend the type of trend. See details.
#' @param alpha the smoothing parameter for constant component. The default is \code{0.2}.
#' @param beta the smoothing parameter for linear component. The default is \code{0.10557}.
#' @param gamma the smoothing parameter for quadratic component. The default is \code{0.07168}.
#' @param lead the number of steps ahead for which prediction is required. 
#' The default is \code{0}.
#' @param plot a logical value indicating to print the plot of original data v.s smoothed 
#' data. The default is \code{TRUE}.
#' @details Simple exponential smoothing is a weighted average between the most recent 
#' observation and the most recent forecasting, with weights \eqn{\alpha} and 
#' \eqn{1 - \alpha}, respectively. To be precise, the smoothing equation of single exponential 
#' smoothing (constant model, \code{trend = 1}) is given by 
#' \deqn{level[t] = \alpha *x[t] + (1 - \alpha)*level[t-1],}
#' and the forecasting equation is 
#' \deqn{hat{x}[t+1|t] = level[t],}
#' for \eqn{t = 1,...,n}.
#' The initial value \eqn{level[0] = x[1]}. For example, \eqn{hat{x}[1|0] = level[0]},
#' \eqn{hat{x}[2|1] = level[1]},..., etc. 
#' 
#' Let \eqn{x1[t]} be the smoothed values of single exponential smoothing. The double
#' exponential smoothing (\code{trend = 2}, a linear model) is to apply a single 
#' exponential smoothing again to the smoothed sequence \eqn{x1[t]}, with a new smoothing 
#' parameter \code{beta}. Similarly, we denote the smoothed values of double 
#' exponential smoothing to be \eqn{x2[t]}. The triple exponential smoothing 
#' (\code{trend = 3}, a quadratic model) is to apply the single exponential smoothing
#' to the smoothed sequence \eqn{x2[t]} with a new smoothing parameter \code{gamma}. The 
#' default smoothing parameters (weights) \code{alpha}, \code{beta}, \code{gamma} are 
#' taken from the equation \code{1 - 0.8^{1/trend}} respectively, which is similar 
#' to the FORECAST procedure in SAS.  
#' 
#' @note Missing values are removed before the analysis.
#' 
#' @return A list with class \code{"es"} containing the following components:
#' \item{estimate}{the smoothed values.}
#' \item{pred}{the predicted values when \code{lead} > 0.}
#' \item{accurate}{the accurate measurements.}
#' @author Debin Qiu
#' @seealso \code{\link{Winters}}, \code{\link{Holt}}, \code{\link{MA}}
#' @examples x <- rnorm(100)
#' es <- expsmooth(x) # trend = 1: a constant model
#' plot(x,type = "l")
#' lines(es$estimate,col = 2)
#' expsmooth(x,trend = 2) # trend = 2: a linear model
#' expsmooth(x,trend = 3) # trend = 3: a quadratic model
#' @importFrom stats ts
#' @importFrom stats frequency
#' @importFrom stats is.ts
#' @export
expsmooth <- function(x,trend = 1,alpha = 0.2,beta = 0.10557,
                        gamma = 0.07168,lead = 0,plot = TRUE)
{
  if (NCOL(x) > 1L)
    stop("'x' must be a univariate time series")
  if (any(c(alpha,beta,gamma) > 1) || any(c(alpha,beta,gamma) < 0))
    stop("'alpha' or 'beta' or 'gamma' must be between 0 and 1")
  if (is.ts(x))
    x <- ts(x[is.finite(x)],start = time(x)[1],frequency = frequency(x))
  else
    x <- x[is.finite(x)]
  if (length(x) < 1L)
    stop("'x' must have at least one observation")
  single.expo <- function(x,parameter) {
    m <- length(x)
    x.hat <- c(x[1],numeric(m - 1))
    for (i in 2:m)  
      x.hat[i] <- parameter*x[i-1] + (1 - parameter)*x.hat[i - 1]
    return(x.hat)
  }
  if (trend == 1) 
    x.hat <- single.expo(x,alpha)
  else if (trend == 2) 
    x.hat <- single.expo(single.expo(x,alpha),beta)
  else if (trend == 3) 
    x.hat <- single.expo(single.expo(single.expo(x,alpha),beta),gamma)
  else stop("'trend' must be one of 1,2,3")
  if (is.ts(x))
    x.hat <- ts(x.hat,start = time(x),frequency = frequency(x))
  result <- list(estimate = x.hat)
  if (lead > 0) {
    if (lead%%1 != 0)
      stop("'lead' must be a positive number")
    result <- c(result,list(pred = rep(x.hat[length(x)],lead)))
  }
  if (plot) {
    plot(x,main = "original v.s smoothed data",type = "l")
    lines(x.hat,col = 2)
  }
  result <- c(result,list(accurate = accurate(x,x.hat,trend,output = FALSE)))
  class(result) <- "es"
  return(result)
}