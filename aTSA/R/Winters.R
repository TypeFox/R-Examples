#' Winters Three-parameter Smoothing
#' @description Performs Winters three-parameter smoothing for a univariate time series
#' with seasonal pattern.
#' @param x a univariate time series.
#' @param period seasonal period. The default is \code{NULL}.
#' @param trend the type of trend component, can be one of 1,2,3 which represents constant,
#' linear and quadratic trend, respectively. The default is \code{trend = 2}.
#' @param lead the number of steps ahead for which prediction is required. 
#' The default is \code{0}.
#' @param plot a logical value indicating to display the smoothed graph. The default is 
#' \code{TRUE}.
#' @param seasonal character string to select an "\code{additive}" or "\code{multiplicative}"
#' seasonal model. The default is "\code{additive}".
#' @param damped a logical value indicating to include the damped trend, only valid for
#' \code{trend = 2}. The default is \code{FALSE}.
#' @param alpha the parameter of level smoothing. The default is \code{0.2}.
#' @param beta the parameter of trend smoothing. The default is \code{0.1057}.
#' @param gamma the parameter of season smoothing. The default is \code{0.07168}.
#' @param phi the parameter of damped trend smoothing, only valid for \code{damped = TRUE}.
#' The default is \code{0.98}.
#' @param a.start the starting value for level smoothing. The default is \code{NA}.
#' @param b.start the starting value for trend smoothing. The default is \code{NA}.
#' @param c.start the starting value for season smoothing. The default is \code{NA}.
#' 
#' @details The Winters filter is used to decompose the trend and seasonal components by 
#' updating equations. This is similar to the function \code{\link{HoltWinters}} in 
#' \code{stats} package but may be in different perspective. To be precise, it uses the 
#' updating equations similar to exponential
#' smoothing to fit the parameters for the following models when 
#' \code{seasonal = "additive"}.
#' If the trend is constant (\code{trend = 1}):
#' \deqn{x[t] = a[t] + s[t] + e[t].}
#' If the trend is linear (\code{trend = 2}):
#' \deqn{x[t] = (a[t] + b[t]*t) + s[t] + e[t].}
#' If the trend is quadratic (\code{trend = 3}):
#' \deqn{x[t] = (a[t] + b[t]*t + c[t]*t^2) + s[t] + e[t].}
#' Here \eqn{a[t],b[t],c[t]} are the trend parameters, \eqn{s[t]} is the seasonal 
#' parameter for the
#' season corresponding to time \eqn{t}.
#' For the multiplicative season, the models are as follows.
#' If the trend is constant (\code{trend = 1}):
#' \deqn{x[t] = a[t] * s[t] + e[t].}
#' If the trend is linear (\code{trend = 2}):
#' \deqn{x[t] = (a[t] + b[t]*t) * s[t] + e[t].}
#' If the trend is quadratic (\code{trend = 3}):
#' \deqn{x[t] = (a[t] + b[t]*t + c[t]*t^2) * s[t] + e[t].}
#' When \code{seasonal = "multiplicative"}, the updating equations for each parameter can 
#' be seen in page 606-607 of PROC FORECAST document of SAS. Similarly, for the 
#' additive seasonal model, the 'division' (/) for \eqn{a[t]} and \eqn{s[t]} in page 606-607 
#' is changed to 'minus' (-).
#' 
#' The default starting values for \eqn{a,b,c} are computed by a time-trend regression over
#' the first cycle of time series. The default starting values for the seasonal factors are 
#' computed from seasonal averages. The default smoothing parameters (weights) \code{alpha, 
#' beta, gamma} are taken from the equation \code{1 - 0.8^{1/trend}} respectively. You can 
#' also use the \code{\link{HoltWinters}} function to get the optimal smoothing parameters
#' and plug them back in this function.
#' 
#' The prediction equation is \eqn{x[t+h] = (a[t] + b[t]*t)*s[t+h]} for \code{trend = 2} and
#' \code{seasonal = "additive"}. Similar equations can be derived for the other options. When
#' the \code{damped = TRUE}, the prediction equation is 
#' \eqn{x[t+h] = (a[t] + (\phi + ... + \phi^(h))*b[t]*t)*s[t+h]}. More details can be 
#' referred to R. J. Hyndman and G. Athanasopoulos (2013).
#' @note The sum of seasonal factors is equal to the seasonal period. This restriction is
#' to ensure the identifiability of seasonality in the models above.
#' @return  A list with class "\code{Winters}" containing the following components:
#' \item{season}{the seasonal factors.}
#' \item{estimate}{the smoothed values.}
#' \item{pred}{the prediction, only available with \code{lead} > 0.}
#' \item{accurate}{the accurate measurements.}
#' @author Debin Qiu
#' @references 
#' P. R. Winters (1960) Forecasting sales by exponentially weighted moving averages, 
#' \emph{Management Science} 6, 324-342.
#' 
#' R. J. Hyndman and G. Athanasopoulos, "Forecasting: principles and practice," 2013. 
#' [Online]. Available: \url{http://otexts.org/fpp/}.
#' @seealso \code{\link{HoltWinters}}, \code{\link{Holt}}, \code{\link{expsmooth}}
#' @examples 
#' Winters(co2)
#' 
#' Winters(AirPassengers, seasonal = "multiplicative")
#' 
#' @importFrom stats coef
#' @importFrom stats lm
#' @importFrom stats ts
#' @importFrom stats plot.ts
#' @importFrom stats frequency
#' @importFrom stats time
#' @export
Winters <- function(x, period = NULL, trend = 2, lead = 0, plot = TRUE,
                    seasonal = c("additive","multiplicative"), damped = FALSE,
                    alpha = 0.2, beta = 0.1057, gamma = 0.07168, phi = 0.98,
                    a.start = NA,b.start = NA,c.start = NA)
{
  DNAME <- deparse(substitute(x))
  if (NCOL(x) > 1L) 
    stop("'x' must be a univariate time series")
  if (any(c(alpha,beta,gamma,phi) > 1) || any(c(alpha,beta,gamma,phi) < 0))
    stop("'alpha', 'beta', 'gamma' or 'phi' must be between 0 and 1")
  if (any(!is.finite(x))) 
    stop("missing values in 'x'")
  if (damped && trend != 2)
    stop("damped trend is only for linear trend component")
  if (is.null(period)) period <- frequency(x) 
  if (!is.ts(x)) x <- as.ts(x[is.finite(x)])
  seasonal <- match.arg(seasonal)
  operat <- function(a,b) 
    switch(seasonal, additive = a - b, multiplicative = a/b)
  n <- length(x)
  s.index <- rep_len(1:period,length.out = n)
  st <- numeric(period)
  for (i in 1:period)  
    st[i] <- mean(x[s.index == i])/mean(x)
  names(st) <- 1:period
  result <- list(season = st)
  phi <- ifelse(damped,phi,1)
  at <- c(ifelse(is.na(a.start),operat(x[1],st[1]),a.start),numeric(n - 1))
  x.hat <- c(x[1],numeric(n - 1))
  s.ind <- function(a) ifelse(a%%period > 0,a%%period,period)
  s.col <- function(b) (0:floor(b/period))*period + b%%period
  if (trend == 1) {   
    for (j in 2:n) { 
      at[j] <- alpha*operat(x[j],st[s.ind(j)]) + (1 - alpha)*at[j - 1]
      st[s.ind(j)] <- gamma*mean(operat(x[s.col(j)],at[s.col(j)])) + (1 - gamma)*st[s.ind(j)]
      x.hat[j] <- switch(seasonal, additive = at[j] + st[s.ind(j)],
                                   multiplicative = at[j]*st[s.ind(j)])
    }
  }
  else if (trend == 2) {
    bt <- c(ifelse(is.na(b.start),operat(x[2],st[2]) - operat(x[1],st[1]),b.start),
            numeric(n - 1))
    for (j in 2:n) {
      at[j] <- alpha*operat(x[j],st[s.ind(j)]) + (1 - alpha)*(at[j - 1] + phi*bt[j - 1])
      bt[j] <- beta*(at[j] - at[j - 1]) + (1 - beta)*phi*bt[j - 1]
      st[s.ind(j)] <- gamma*mean(operat(x[s.col(j)],at[s.col(j)] + phi*bt[s.col(j)])) + 
                      (1 - gamma)*st[s.ind(j)]
      x.hat[j] <- switch(seasonal, additive = at[j] + phi*bt[j] + st[s.ind(j)],
                         multiplicative = (at[j] + bt[j])*st[s.ind(j)])
    }
  }
  else if (trend == 3) {
    tt <- 1:period
    yy <- x[1:period]/st
    coef0 <- coef(lm(yy ~ tt + I(tt^2)))
    bt <- c(ifelse(is.na(b.start),coef0[2],b.start),numeric(n - 1))
    ct <- c(ifelse(is.na(c.start),coef0[3],c.start),numeric(n - 1))
    for (j in 2:n) {
      at[j] <- alpha*operat(x[j],st[s.ind(j)]) + 
               (1 - alpha)*(at[j - 1] + bt[j - 1] + ct[j - 1])
      bt[j] <- beta*(at[j] - at[j - 1] + ct[j - 1]) + 
               (1 - beta)*(bt[j - 1] + 2*ct[j - 1])
      ct[j] <- beta*(bt[j] - bt[j - 1])/2 + (1 - beta)*ct[j - 1]
      st[s.ind(j)] <- gamma*mean(operat(x[s.col(j)],at[s.col(j)] + bt[s.col(j)] + 
                                       ct[s.col(j)])) + (1 - gamma)*st[s.ind(j)]
      x.hat[j] <- switch(seasonal, additive = at[j] + bt[j] + ct[j] + st[s.ind(j)],
                         multiplicative = (at[j] + bt[j] + ct[j])*st[s.ind(j)])
    }
  }
  else stop("'trend' must be one of 1,2,3")
  x.hat <- ts(x.hat,start = time(x)[1],frequency = frequency(x))
  result <- c(result,list(estimate = x.hat))
  if (lead > 0) {
    index <- (n + 1:lead)%%period
    if (trend == 1) 
      x.pred <- switch(seasonal,additive = rep(at[n],lead) + st[index],
                                multiplicative = rep(at[n],lead) * st[index])
    if (trend == 2) 
      x.pred <- switch(seasonal,additive = rep(at[n],lead) + cumsum(phi^(1:lead))*
                         rep(bt[n],lead)*(n + 1:lead)/n + st[index],
                       multiplicative = (rep(at[n],lead) + cumsum(phi^(1:lead))*
                                           rep(bt[n],lead)*(n + 1:lead)/n) * st[index])
    if (trend == 3) 
      x.pred <- switch(seasonal,additive = rep(at[n],lead) + rep(bt[n],lead)*
                        (n + 1:lead)/n + rep(ct[n],lead)*((n + 1:lead)/n)^2 + st[index],
                        multiplicative = (rep(at[n],lead) + rep(bt[n],lead)*
                        (n + 1:lead)/n + rep(ct[n],lead)*((n + 1:lead)/n)^2) * st[index])
    pred <- matrix(c((n + 1:lead),(n + 1:lead)%%period,x.pred),nrow = 3,byrow = TRUE)
    dimnames(pred) <- list(c("time:","period:","pred:"),c(rep("",lead)))
    result <- c(result,list(pred = pred))
  }
  if (plot) {
    plot.ts(x,main = "original v.s smoothed data",ylab = DNAME)
    lines(x.hat,col = 2)
  }
  result <- c(result,list(accurate = accurate(x,x.hat,trend,output = FALSE)))
  class(result) <- "Winters"
  return(result)
}