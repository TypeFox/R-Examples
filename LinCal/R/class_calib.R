#' Classical Linear Calibration Function
#'
#' \code{class.calib} uses the classical frequentist approach to estimate an unknown X given observed vector y0 and calculates confidence interval estimates.
#' @param x numerical vector of regressor measurments
#' @param y numerical vector of observation measurements
#' @param alpha the confidence interval to be calculated
#' @param y0 vector of observed calibration value
#' @references Eisenhart, C. (1939). The interpretation of certain regression methods and their use in biological and industrial research. Annals of Mathematical Statistics. 10, 162-186.
#' @keywords linear calibration
#' @export
#' @examples
#' X <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
#' Y <- c(1.8,1.6,3.1,2.6,3.6,3.4,4.9,4.2,6.0,5.9,6.8,6.9,8.2,7.3,8.8,8.5,9.5,9.5,10.6,10.6)
#'
#' class.calib(X,Y,0.05,6)

class.calib <- function(x, y, alpha, y0){



  reg1   <- lm(y ~ x)
  x.pre  <- (mean(y0) - coef(reg1)[1])/coef(reg1)[2]



  alpha1 = 2*alpha
  n      <- length(y)
  x.mean <- mean(x)

  Sxx   <- sum(x^2) - (sum(x)^2)/n
  s2.e  <- sum(reg1$resid^2)/(n - 2)
  s.e   <- sqrt(s2.e)
  tval  <- qt(1-alpha1, df = n-2, ncp=0, lower.tail = TRUE, log.p = FALSE)
  t.obs <- coef(reg1)[2]/(s.e/sqrt(Sxx))
  c2    <-  round(((tval^2)*s2.e) / (coef(reg1)[2]^2 * Sxx),4)
  d     <- ((tval * s.e)/coef(reg1)[2]) * sqrt(((n+1)/n) * (1 - c2) + ((x.pre - x.mean)^2)/Sxx)

  upper <- x.mean + (1/(1 - c2))*((x.pre - x.mean) + d)
  lower <- x.mean + (1/(1 - c2))*((x.pre - x.mean) - d)
  limits   <- cbind(min(lower,upper),max(lower,upper))

  list(x.pre = round(x.pre,9), lim=limits)

}
