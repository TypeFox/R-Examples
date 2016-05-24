#' Bayesian Classical Linear Calibration Function
#'
#' \code{huntlam.calib} uses the classical Bayesian approach to estimate an unknown X given observed vector y0 and calculates credible interval estimates.
#' @param x numerical vector of regressor measurments
#' @param y numerical vector of observation measurements
#' @param alpha the confidence interval to be calculated
#' @param y0 vector of observed calibration value
#' @references Hunter, W., and Lamboy, W. (1981). A Bayesian Analysis of the Linear Calibration Problem. Technometrics. 3, 323-328.
#' @keywords linear calibration
#' @export
#' @examples
#' X <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
#' Y <- c(1.8,1.6,3.1,2.6,3.6,3.4,4.9,4.2,6.0,5.9,6.8,6.9,8.2,7.3,8.8,8.5,9.5,9.5,10.6,10.6)
#'
#' huntlam.calib(X,Y,0.05,6)

huntlam.calib <- function(x,y, alpha, y0){


  reg4  <- lm(y ~ x)
  x.hl  <- (mean(y0) - reg4$coef[1])/reg4$coef[2]

  anova4 <- anova(reg4)
  n      <- length(y)
  s2     <- anova4$Mean[2]

  des1   <- cbind(1,x)
  S.mat  <- s2 * solve(t(des1) %*% des1)
  V.mat  <- matrix(c(S.mat[1,1]+(s2),-S.mat[1,2], -S.mat[1,2], S.mat[2,2]),2,2)
  scale  <- (V.mat[1,1]*V.mat[2,2]-V.mat[1,2]**2)/(V.mat[2,2]*reg4$coef[2]**2)

  alpha4 = 1*alpha
  x.pre   <- matrix(0,length(mean(y0)))
  quant   <- matrix(0,length(mean(y0)),2)
  for(k in 1:length(mean(y0))){
    x.pre[k]    <- rnorm(1,x.hl[k],sqrt(scale/100))
    quant[k,]   <- qnorm(c(alpha4, (1-alpha4)),mean = x.hl[k], sd=sqrt(scale))
  }

  list(x.pre=x.hl, sd=sqrt(scale), lim=quant)

}

