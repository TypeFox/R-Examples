#' Identify a Time Series Model
#' @description Checks the white noise and stationarity of a univariate time series, 
#' and also identifies an appropriate ARIMA model using AICC criterion.
#' @param x a numeric vector or univariate time series.
#' @param p the maximum lag order for AR process. The default is \code{NULL}.
#' @param q the maximum lag order for MA process. The default is \code{NULL}.
#' @param nlag the lag parameter to calculate the Ljung-Box test statistic. 
#' The default is \code{6}.
#' @param intercept an intercept to be included in ARIMA model, only valid for 
#' \code{p} > 0 or \code{q} > 0. The default is \code{TRUE}.
#' @param stat.test the stationary test for time series, see \code{\link{stationary.test}} for
#' more details. The default is \code{FALSE}.
#' @param method the method of stationary test, only valid for \code{stat.test = TRUE}.
#' @param output a logical value indicating to print the results in R console. 
#' The default is \code{TRUE}.
#' 
#' @details This function is similar to IDENTIFY statement in ARIMA procedure of SAS software, 
#' which is to check the white noise and stationarity for a univariate time 
#' series. The white noise check is accomplished by using \code{\link{Box.test}} 
#' in \code{stats} package, with the default method \code{type = "Ljung-Box"}. 
#' The stationary check uses the \code{\link{stationary.test}} implemented in this package. 
#' 
#' The AICC criterion (Burnham and Anderson (2002)) is used to identify an optimal 
#' model which has the minimum AICC value. The AICC is defined as 
#' \deqn{AICC = AIC + 2k(k+1)/(n-k-1),}
#' where \eqn{AIC = 2k - 2log(Loglik)} which is called Akaike information criterion 
#' (Akaike (1974)). Here, 
#' \eqn{k,n} are the number of estimated parameters and observations, respectively. 
#' \eqn{Loglik} is the maximized value of the likelihood function for the model.
#' 
#' Four plots are made: plot of original data, ACF plot, PACF plot and p.value of white 
#' noise check. 
#' @note Missing values are removed before the analysis.
#'  
#' @return A list with class "\code{identify}" containing the following components:
#'  \item{WNcheck}{a matrix with three columns for results of white noise check for each lag.}
#'  \item{aicc}{the minimum AICC value, only available for \code{p > 0} or \code{q > 0}.}
#'  \item{min.p}{the optimal order \code{p} for AR process, only available for 
#'  \code{p > 0} or \code{q > 0}.}
#'  \item{min.q}{the optimal order \code{q} for AR process, only available for 
#'  \code{p > 0} or \code{q > 0}.}
#'  \item{stnt.test}{a list of stationary test results with three components. See 
#'  \code{\link{stationary.test}} for more details.}
#'  
#' @author Debin Qiu
#' @references
#' Akaike, H. (1974), "A new look at the statistical model identification", \emph{IEEE 
#' Transactions on Automatic Control}, 19 (6): 716-723.
#' 
#' Box, G. E. P. and Pierce, D. A. (1970), Distribution of residual correlations 
#' in autoregressive-integrated moving average time series models. \emph{Journal 
#' of the American Statistical Association}, 65, 1509-1526.
#' 
#' Burnham, K. P.; Anderson, D. R. (2002), Model Selection and Multimodel Inference: 
#' A Practical Information-Theoretic Approach (2nd ed.), Springer-Verlag
#'  
#' Ljung, G. M. and Box, G. E. P. (1978). On a measure of lack of fit in time series 
#' models. \emph{Biometrika} 65, 297-303.
#' 
#' Harvey, A. C. (1993) Time Series Models. 2nd Edition, Harvester Wheatsheaf, 
#' NY, pp. 44, 45.
#' 
#' @examples x <- arima.sim(list(order = c(2,0,0),ar = c(0.2,0.4)),n = 100)
#' identify(x) # white noise check
#' identify(x,stat.test = TRUE) # white noise and stationarity check
#' identify(x,p = 3,q = 2) # white noise check and optimal model identification.
#' @importFrom stats Box.test
#' @importFrom stats acf
#' @importFrom stats pacf
#' @importFrom graphics abline
#' @importFrom stats plot.ts
#' @importFrom stats arima
#' @importFrom graphics par
#' @importFrom graphics plot
#' @importFrom stats sd
#' @export 
identify <- function(x,p = NULL,q = NULL,nlag = 6,intercept = TRUE, 
                     stat.test = FALSE,method = c("adf","pp","kpss"),
                     output = TRUE)
{
  DNAME <- deparse(substitute(x))
  if (NCOL(x) > 1)
    stop("'x' must be a numeric vector or univariate time series")
  x <- x[is.finite(x)]
  n <- length(x)
  if (n < 1L)
    stop("invalid length of 'x'")
  mu <- mean(x)
  SD <- sd(x)
  QB <- p.value <- NULL
  for (i in 1:nlag) {
    LB.test <- Box.test(x,lag = i, type = "Ljung-Box")
    QB[i] <- LB.test$statistic
    p.value[i] <- LB.test$p.value
  }
  WNcheck <- cbind(1:nlag,QB,p.value)
  colnames(WNcheck) <- c("lag","LB","p.value")
  result <- list(WNcheck = WNcheck)
  op <- par(mfrow = c(2,2))
  plot.ts(x,main = paste("plot of",DNAME))
  acf(x,lag.max = nlag,main = paste("ACF plot of", DNAME))
  pacf(x,lag.max = nlag,main = paste("PACF plot of", DNAME))
  plot(1:nlag,p.value,xlab = "Lag", ylim = c(0,1),
       main = "p.value of WN check",pch = 19,col = 2)
  abline(0.05,0,lty = 2, col = 3)
  par(op) 
  if (!is.null(p) || !is.null(q)) {
    p <- ifelse(is.null(p),0,p)
    q <- ifelse(is.null(q),0,q)
    if ((p < 0 || q < 0) || (p%%1 != 0 || q%%1 != 0))
      stop("'p' or 'q' must be a positive integer")
    p.s <- 0:p
    q.s <- 0:q
    aicc <- matrix(0,p + 1,q + 1)
    for (i in p.s) {
      for (j in q.s) {
        aic <- arima(x,order = c(i,0,j),include.mean = intercept)$aic
        k <- max(i + j + ifelse(intercept,1,0),1)
        aicc[i + 1,j + 1] <-  aic + 2*k*(k + 1)/(n - k - 1)
      }      
    }
    rownames(aicc) <- paste("AR",p.s,sep = "")
    colnames(aicc) <- paste("MA",q.s,sep = "")
    index <- which(aicc == min(aicc), arr.ind = TRUE)
    result <- c(result,list(aicc = min(aicc), min.p = index[1] - 1,
                   min.q = index[2] - 1))
  }
  if (output) {
    cat("Name of Variable:",DNAME,"\n")
    cat("Number of Observations:", n, "\n")
    cat("Mean of Working Series:",mu, "\n")
    cat("Standard Deviation:", SD, "\n")
    cat("--------------------","\n")
    cat("Autocorrelation Check for White Noise","\n")
    cat("alternative: non-white-noise","\n\n")
    print(WNcheck,digits = 3) 
    if (!is.null(p) || !is.null(q)) {
      cat("--------------------","\n")
      cat("Minimum Information Criterion:","\n")  
      print(aicc,digits = 5)
      cat("\n")
      cat("Minimum Table Value:",paste("AICC(",index[1] - 1,",",
           index[2] - 1,")", sep = "")," = ", aicc[index],"\n")
    }
  }
  if (stat.test) {
    cat("---------------------","\n")
    method <- match.arg(method)
    METHOD <- switch(method,adf = "adf",pp = "pp",kpss = "kpss")
    result <- c(result,list(stn.test = stationary.test(x,method = METHOD,
                                                       output = output)))
  }
  class(result) <- "identify"
  identify <- result
}