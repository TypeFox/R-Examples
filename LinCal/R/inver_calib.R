#' Inverse Linear Calibration Function
#'
#' \code{inver.calib} uses the inverse frequentist approach to estimate an unknown X given observed vector y0 and calculates confidence interval estimates.
#' @param x numerical vector of regressor measurments
#' @param y numerical vector of observation measurements
#' @param alpha the confidence interval to be calculated
#' @param y0 vector of observed calibration value
#' @references Krutchkoff, R. G. (1967). Classical and Inverse Regression Methods of Calibration. Technometrics. 9, 425-439.
#' @keywords linear calibration
#' @export
#' @examples
#' X <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
#' Y <- c(1.8,1.6,3.1,2.6,3.6,3.4,4.9,4.2,6.0,5.9,6.8,6.9,8.2,7.3,8.8,8.5,9.5,9.5,10.6,10.6)
#'
#' inver.calib(X,Y,0.05,6)


inver.calib <- function(x,y, alpha, y0){



  reg2   <- lm(x~y)
  x.pre2 <- coef(reg2)[1] + coef(reg2)[2] * mean(y0)


  alpha2 = 2*alpha
  anova1 <- anova(reg2)
  n      <- length(y)
  s      <- anova1$Mean[2]

  confint <- matrix(0, length(mean(y0)),2)
  for(j in 1:length(mean(y0))){
    confint[j,]  <- qnorm(c(alpha2,(1-alpha2)), x.pre2[j], sqrt(s))
  }


  list(x.pre=round(x.pre2,9), lim=confint)

}

