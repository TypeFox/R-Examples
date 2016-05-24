#' Stepwise Autoregressive Model
#' @description Fit a stepwise autoregressive model
#' @param y a numeric vector of response
#' @param xreg a numeric vector or matrix of exogenous input variables. The default is 
#' \code{NULL}.
#' @param trend the type of trend with respective to time. The default is \code{linear}.
#' @param order the order to fit the AR model for residuals. The default is \code{NULL}.
#' @param lead the number of steps ahead for which prediction is required. 
#' The default is \code{0}.
#' @param newx a matrix of new data of \code{xreg} for predictions. The default is
#' \code{NULL}.
#' @param output a logical value indicating to print the results in R console. The default is 
#' \code{NULL}.
#' @param ... additional arguments for \code{\link{ar}} function.
#' 
#' @details The stewise autoregressive model uses a two-stage procedure to fit time series.
#' The first stage is to fit a (\code{constant},\code{linear},\code{quadratic})
#' model with respective to time sequence: 
#' \eqn{t = (1:n)/n}, where \eqn{n = length(y)}. If \code{xreg} is supplied, 
#' the fitted model is updated by 
#' \deqn{y = \mu + \beta*xreg + e[t] }
#' for \code{trend = "constant"}, and 
#' \deqn{y = \mu + \beta*xreg + \alpha*t + e[t]}
#' for \code{trend = "linear"}, and 
#' \deqn{y = \mu + \beta*xreg + \alpha[1]*t + \alpha[2]*t^2 + e[t]}
#' for \code{trend = "quadratic"}.
#' The second stage is to fit an autoregressive process to the residuals of the fitted
#' model obtained in the first stage, which is accomplished by using \code{\link{ar}} function 
#' in \code{stats} package.
#' 
#' @note If \code{lead} > 0 and \code{xreg} is supplied, \code{newx} must also be 
#' supplied in order to make a prediction. The \code{newx} must be a matrix with 
#' the same number of columns as \code{xreg} and the number of rows being 
#' equal to \code{lead}. The predictions should be used with cautions.
#' @return A list with class "\code{stepar}" containing the following components:
#' \item{coef}{a estimated coefficient matrix including the t-test results.}
#' \item{sigma}{the square root of the estimated variance of the random error.}
#' \item{R.squared}{the R^2 for fitted model in the first stage.}
#' \item{pred}{the predictions, only available for \code{lead} > 0.}
#' @author Debin Qiu
#' @examples x <- 5*(1:100)/100
#' x <- x + arima.sim(list(order = c(1,0,0),ar = 0.4),n = 100)
#' stepar(x)
#' stepar(x,order = 1)
#' 
#' # with xreg supplied
#' X <- matrix(rnorm(200),100,2)
#' y <- 0.1*X[,1] + 1.2*X[,2] + rnorm(100)
#' stepar(y,X)
#' # make a prediction with lead = 1; used with caution.
#' newdat1 <- matrix(rnorm(2),nrow = 1)
#' fit1 <- stepar(y,X,lead = 1,newx = newdat1,output = FALSE)
#' # make a prediction with lead = 2; used with caution.
#' newdat2 <- matrix(rnorm(4),nrow = 2)
#' fit2 <- stepar(y,X,lead = 2,newx = newdat2,output = FALSE)
#' @importFrom stats lm
#' @importFrom stats ar
#' @importFrom stats residuals
#' @importFrom stats pt
#' @importFrom stats predict
#' @importFrom stats coef
#' @export
stepar <- function(y,xreg = NULL, trend = c("linear","quadratic","constant"),
                   order = NULL, lead = 0, newx = NULL, output = TRUE,...)
{
  if (NCOL(y) > 1L)
    stop("'y' must be a numeric vector or univariate time series")
  if (lead < 0 || lead%%1 != 0)
    stop("'lead' must be a positive integer")
  if (!is.null(xreg) && NROW(y) != NROW(xreg)) 
    stop("'y' and 'xreg' must have same length")  
  trend <- match.arg(trend)
  xreg <- xreg[is.finite(y),]
  y <- y[is.finite(y)]
  n <- length(y)
  if (n < 1L)
    stop("invalid length of 'y'")
  t <- (1:n)/n
  if (!is.null(xreg)) 
    model <- switch(trend,constant = lm(y ~ xreg) , 
                    linear = lm(y ~ xreg + t),
                    quadratic = lm (y ~ xreg + t + I(t^2)))
  else 
    model <- switch(trend,constant = lm(y ~ 1), linear = lm(y ~ t), 
                    quadratic = lm(y ~ t + I(t^2)))
  res <- residuals(model)
  coef.mtrx <- as.matrix(summary(model)$coefficients)
  ar.fit <- if (!is.null(order)) 
                ar(res,aic = FALSE, order.max = order,...) 
            else ar(res,aic = TRUE, order.max = NULL,...)
  order <- ar.fit$order
  if (order > 0) {
    coef.ar <- ar.fit$ar
    se.ar <- sqrt(diag(ar.fit$asy.var.coef))
    t.value <- coef.ar/se.ar
    p.value <- pt(t.value,n - length(coef.ar))
    p.value <- 2*pmin(p.value, 1 - p.value)
    ar.coef.mtrx <- matrix(c(coef.ar,se.ar,t.value,p.value),ncol = 4)
    rownames(ar.coef.mtrx) <- paste("AR",1:order,sep = "")
    coef.mtrx <- rbind(coef.mtrx,ar.coef.mtrx) 
  }
  sigma2 <- sum(res^2)/(n - 1)
  R2 <- summary(model)$r.squared
  result <- list(coef = coef.mtrx,sigma = sqrt(sigma2),R.squared = R2)
  if (lead > 0) {
    newdata <- switch(trend,constant = rep(1,lead),
                      linear = cbind(rep(1,lead),(n + 1:lead)/n),
                      quadratic = cbind(rep(1,lead),(n + 1:lead)/n,
                                        ((n + 1:lead)/n)^2))
    newdata <- matrix(newdata,nrow = lead)
    if (!is.null(xreg)) {
      if (is.null(newx))
        stop("please provide 'newx' for prediction")
      else {
        if (!is.matrix(newx))
          stop("'newx' must be a supplied matrix")
        else {
          if (NCOL(newx) != NCOL(xreg))
            stop(paste("NCOL(newx) must be equal to",NCOL(xreg)))
          if (NROW(newx) != lead)
            stop(paste("NROW(newx) must be equal to",lead))
        }
      }
      newdata <- matrix(cbind(newdata, newx),nrow = lead)
    }
    pred <- newdata%*%matrix(coef(model),ncol = 1)
    if (order > 0) {
      pred.ar <- predict(ar.fit,n.ahead = lead)
      pred <- pred + pred.ar$pred    
    }
    result <- c(result,list(pred = pred))
  }
  if (output) {
    cat("Parameter of estimates for stepwise AR model:","\n")
    print(coef.mtrx,digits = 3)
    cat("------","\n")
    cat("sigma^2 estimated as:",sigma2,";","R.squared =",R2)
  }
  class(result) <- "stepar"
  stepar <- result
}