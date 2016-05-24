#' ARCH Engle's Test for Residual Heteroscedasticity
#' @description Performs Portmanteau Q and Lagrange Multiplier tests for the null 
#' hypothesis that the residuals of a ARIMA model are homoscedastic.
#' @param object an object from arima model estimated by 
#' \code{\link{arima}} or \code{\link{estimate}} function.
#' @param output a logical value indicating to print the results in R console, including the
#' plots. The default is \code{TRUE}.
#' @details The ARCH Engle's test is constructed based on the fact that if the residuals 
#' (defined as \eqn{e[t]}) are heteroscedastic, the squared residuals (\eqn{e^2[t]}) are 
#' autocorrelated. The 
#' first type of test is to examine whether the squares of residuals are a sequence of white 
#' noise, which is called Portmanteau Q test and similar to the Ljung-Box test on the squared 
#' residuals. The second type of test proposed by Engle (1982) is the Lagrange Multiplier 
#' test which is to fit a linear regression model for the squared residuals and examine 
#' whether the fitted model
#' is significant. So the null hypothesis is that the squared residuals are a sequence
#' of white noise, namely, the residuals are homoscedastic. The lag parameter
#' to calculate the test statistics is taken from an integer sequence of \eqn{1:min(24,n)} with
#' step 4 if \eqn{n > 25}, otherwise 2, where \eqn{n} is the number of nonmissing observations. 
#' 
#' The plots of residuals, squared residuals, p.values of PQ and LM tests will be drawn if 
#' \code{output = TRUE}.
#' @note Missing values are removed before analysis. 
#' @return A matrix with the following five columns:
#' \item{\code{order}}{the lag parameter to calculate the test statistics.}
#' \item{\code{PQ}}{the Portmanteau Q test statistic.}
#' \item{\code{p.value}}{the p.value for PQ test.}
#' \item{\code{LM}}{the Lagrange Multiplier test statistic.}
#' \item{\code{p.value}}{the p.value for LM test.}
#' @author Debin Qiu
#' @references
#' Engle, Robert F. (1982). Autoregressive Conditional Heteroscedasticity with Estimates 
#' of the Variance of United Kingdom Inflation. \emph{Econometrica}, 50 (4): 987-1007.
#' 
#' McLeod, A. I. and W. K. Li. Diagnostic Checking ARMA Time Series Models Using 
#' Squared-Residual Autocorrelations. \emph{Journal of Time Series Analysis}. 
#' Vol. 4, 1983, pp. 269-273.
#' 
#' @examples x <- rnorm(100)
#' mod <- estimate(x,p = 1) # or mod <- arima(x,order = c(1,0,0))
#' arch.test(mod)
#' @importFrom stats residuals
#' @importFrom stats Box.test
#' @importFrom stats embed
#' @importFrom stats lm
#' @importFrom stats pchisq
#' @importFrom graphics plot
#' @importFrom graphics abline
#' @importFrom graphics par
#' @export
arch.test <- function(object,output = TRUE)
{
  if (class(object) != "Arima" && class(object) != "estimate")
    stop("'object' should be 'Arima' or 'estimate' class estimated 
         from arima() or estimate()")
  class(object) <- "Arima"
  res <- residuals(object)
  res2 <- res^2
  n <- length(res)
  if (n < 2L)
    stop("not enough length of residuals")
  lag.max <- min(24,n)
  lag.seq <- seq(ifelse(lag.max > 4,4,2),lag.max, by = 4)
  PQ <- LM <- p.value1 <- p.value2 <- numeric(length(lag.seq))
  for (i in 1:length(lag.seq)) {
    k <- lag.seq[i]
    LB.test <- Box.test(res2,lag = k)
    PQ[i] <- LB.test$statistic
    p.value1[i] <- LB.test$p.value  
    rt <- embed(res2,k)
    lm.fit <- lm(rt[,1] ~ rt[,-1])
    SSE <- sum((rt[,1] - residuals(lm.fit))^2)
    SST <- sum((rt[,1] - mean(rt[,1]))^2)
    LM[i] <- ((SST - SSE)/k)/(SSE/(n - 2*k - 1))
    p.value2[i] <- 1 - pchisq(LM[i],k - 1)
  }
  PQcheck <- matrix(c(lag.seq,PQ,p.value1),ncol = 3)
  LMcheck <- matrix(c(LM,p.value2),ncol = 2)
  result <- cbind(PQcheck,LMcheck)
  colnames(result) <- c("order","PQ","p.value","LM","p.value")
  if (output) {
    cat("ARCH heteroscedasticity test for residuals","\n")
    cat("alternative: heteroscedastic","\n\n")
    cat("Portmanteau-Q test:","\n")
    print(result[,1:3],digits = 3) 
    cat("Lagrange-Multiplier test:","\n")
    print(result[,c(1,4,5)],digits = 3)
    op <- par(mfrow = c(2,2))
    plot(1:n,res,xlab = "Time",ylab = "resid")
    abline(h = 0,col = 2,lty = 2)
    plot(1:n,res2,xlab = "Time",ylab = "resid^2")
    abline(h = 0,col = 2,lty = 2)
    plot(lag.seq,p.value1,ylim = c(0,1),xlab = "Order",ylab = "PQ prob.")
    abline(h = .05,col = 3, lty = 2)
    plot(lag.seq,p.value2,ylim = c(0,1),xlab = "Order",ylab = "LM prob.")
    abline(h = .05,col = 3, lty = 2)
    par(op)
  }
  arch.test <- result 
}