#' @title Register method for cdfqr object functions
#' @description Register method for cdfqr object functions.
#' @aliases residuals.cdfqr
#' @param object The cdfqr model project  
#' @param type The type of residuals to be extracted: \code{'raw'}, \code{'pearson'},\code{'std.pearson'}, or \code{'deviance'},
#' @param ... currently ignored 
#' @return residuals of a specified type.
#' @method residuals cdfqr
#' @export
#' 
#' @examples
#' data(cdfqrExampleData)
#' fit <- cdfquantreg(crc99 ~ vert | confl, 't2','t2', data = JurorData)
#' 
#' residuals(fit, "pearson")

residuals.cdfqr <- function(object, type = c("raw","pearson", "deviance"), ...) {
  
  type <- match.arg(type)
  
  residuals <- object$residuals  # Obtain the raw residuals
  type <- match.arg(type)
  
  ydata <- object$y  # observed data
  n <- length(ydata)  # number of observations
  fitted <- fitted(object,"full")  # model fitted values
  mu <- fitted(object, type = "mu")  # fitted mu values
  sigma <- fitted(object, type = "sigma")  # fitted sigma values values
  
  # Get the distribution and the liklihood function
  dist <- object$family
  fd <- dist$fd
  sd <- dist$sd
 
  residuals <- ydata - fitted  # Raw residuals
  pearson <- residuals/as.numeric(sqrt(var(fitted)))  # Pearson
 
  
  # - Deviance residuals
  deviance_r <- sign(residuals) * sqrt(2 * abs(-qrLogLik(ydata, mu, sigma, fd, sd) + qrLogLik(fitted, 
    mu, sigma, fd, sd)))
  
    res <- switch(type, raw = {
    residuals
  }, pearson = {
    pearson
  },  deviance = {
    deviance_r
  })
  
  return(res)
} 
