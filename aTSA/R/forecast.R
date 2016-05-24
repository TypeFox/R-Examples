#' Forecast From ARIMA Fits
#' @description Forecasts from models fitted by \code{\link{arima}} or \code{\link{estimate}}
#' function.
#' 
#' @param object the result of an \code{arima} or \code{estimate} fit.
#' @param lead the number of steps ahead for which prediction is required. The default is 
#' \code{1}.
#' @param id the id of the observation which is the time. The default is \code{NULL}.
#' @param alpha the significant level for constructing the confidence interval of prediction. 
#' The default is \code{0.05}. 
#' @param output a logical value indicating to print the results in R console. 
#' The default is \code{TRUE}.
#' 
#' @details This function is originally from \code{\link{predict.Arima}} in \code{stats} package, 
#' but has a nice output
#' including 100*(1 - \eqn{\alpha})\% confidence interval and a prediction plot. It is 
#' similar to FORECAST statement in PROC ARIMA of SAS. 
#' @return A matrix with \code{lead} rows and five columns. Each column represents the number 
#' of steps ahead (\code{Lead}), the predicted values (\code{Forecast}), the standard errors 
#' (\code{S.E}) and the 100*(1 - \eqn{\alpha})\% lower bound (\code{Lower}) and upper bound 
#' (\code{Upper}) of confidence interval. 
#' @author Debin Qiu
#' @seealso \code{\link{predict.Arima}}
#' @examples x <- arima.sim(list(order = c(3,0,0),ar = c(0.2,0.4,-0.15)),n = 100)
#' fit <- estimate(x,p = 3) # same as fit <- arima(x,order = c(3,0,0))
#' forecast(fit,lead = 4)
#' 
#'# forecast with id
#'t <- as.Date("2014-03-25") + 1:100
#'forecast(fit,lead = 4, id = t)
#'@importFrom stats predict
#'@importFrom graphics polygon
#'@importFrom stats qnorm
#'@importFrom graphics plot
#'@importFrom graphics lines
#'@export    
forecast <- function(object,lead = 1,id = NULL,alpha = 0.05,output = TRUE)
{
  if (class(object) != "Arima" && class(object) != "estimate")
    stop("'object' should be 'Arima' or 'estimate' class estimated 
         from arima() or estimate()")
  class(object) <- "Arima"
  pred.arima <- predict(object,n.ahead = lead)
  pred <- pred.arima$pred
  se <- pred.arima$se
  n <- length(object$residuals)
  z <- qnorm(1 - alpha/2)
  if (!is.null(id)) {
    if (class(id) != "Date")
      stop("'id' must be 'Date' class")
    n <- id[n]
  }
  nhead <- n + 1:lead
  result <- matrix(c(1:lead, pred,se,pred - z*se, pred + z*se),lead,5)
  colnames(result) <- c("Lead","Forecast","S.E","Lower","Upper")
  rownames(result) <- as.character(nhead)
  if (output) {
    cat("Forecast for univariate time series:","\n")
    print(result,digits = 3)
    cat("------","\n")
    cat("Note: confidence level =", 100*(1 - alpha),"%","\n")
    x <- c(eval(parse(text = object$series)),result[,2])
    id <- c(if (is.null(id)) 1:n else id,nhead) 
    plot(id,x,xlab = "time",ylab = object$series,type = "l",
         ylim = c(min(x,result[,4]),max(x,result[,5])))   
    polygon(c(nhead,rev(nhead)),c(result[,4],rev(result[,5])), col = "gray")
    lines(nhead,result[,2],lty = 1, col = 2,lwd = 2)
    lines(nhead,result[,4],lty = 2, col = 3,lwd = 2)
    lines(nhead,result[,5],lty = 2, col = 3,lwd = 2)
  } 
  forecast <- result
}