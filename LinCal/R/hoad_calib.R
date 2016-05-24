#' Bayesian Inverse Linear Calibration Function
#'
#' \code{hoad.calib} uses an inverse Bayesian approach to estimate an unknown X given observed vector y0 and calculates credible interval estimates.
#' @param x numerical vector of regressor measurments
#' @param y numerical vector of observation measurements
#' @param alpha the confidence interval to be calculated
#' @param y0 vector of observed calibration value
#' @references Hoadley, B. (1970). A Bayesian look at Inverse Linear Regression. Journal of the American Statistical Association. 65, 356-369.
#' @keywords linear calibration
#' @export
#' @examples
#' X <- c(1,1,2,2,3,3,4,4,5,5,6,6,7,7,8,8,9,9,10,10)
#' Y <- c(1.8,1.6,3.1,2.6,3.6,3.4,4.9,4.2,6.0,5.9,6.8,6.9,8.2,7.3,8.8,8.5,9.5,9.5,10.6,10.6)
#'
#' hoad.calib(X,Y,0.05,6)

hoad.calib <- function(x,y,alpha,y0){

  df <- data.frame(x=x,y=y)
  sxx    <- sum((df$x-mean(df$x))^2)
  n      <- length(df$x)
  a      <- sqrt(n/sxx)
  data_scaled   <- df
  data_scaled$x <- a*((df$x)-mean(df$x))

  y      <- summary(lm(y~x,data_scaled))
  anova3 <- anova(lm(y~x,data_scaled))
  sigma  <- anova3$Mean[2]
  z      <- summary(lm(x~y,data_scaled))
  F      <- n*((y$coeff[2])**2)/((y$sigma**2)/1)
  R      <- F/(F+n-2)


  x_inv  <- z$coeff[1] + (mean(y0) * z$coeff[2])
  var1   <- (n+1+((x_inv**2)/R))/(F+n-2)


  alpha3 = 2*alpha
  t      <- qt((1-alpha3),(n-2))
  lower1 <- x_inv - (t*sqrt(var1))
  upper1 <- x_inv + (t*sqrt(var1))
  x_inv  <- (x_inv/a)  + mean(df$x)
  lower  <- (lower1/a) + mean(df$x)
  upper  <- (upper1/a) + mean(df$x)
  sd1    <- sqrt(var1/(a^2))
  limits3<-cbind(lower,upper)
  x.pre  <- matrix(0,length(mean(y0)))
  for(h in 1:length(mean(y0))){
    x.pre[h] <- rnorm(1,x_inv[h],sd1[h]/sqrt(n))
  }
  list(x.pre=x_inv, sd = sd1, lim=limits3)

}

