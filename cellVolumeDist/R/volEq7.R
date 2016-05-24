`volEq7` <-
function(A=1,r=1,sigma_r=1,t=29,sigma_t=.3*t, V=1) {

  intEq <- function(tau) {

    term2 <- 1/3*(t*sigma_r)^2 + (tau*sigma_r)^2
    
    (2^(-tau/t))*(1/(sqrt(1/3*(r*sigma_t)^2 + term2) * sqrt(2*pi))) * exp( -(V - r*(t + tau))^2 / (2* ((1/3*(r*sigma_t)^2 + term2))))* (.5-(.5*erf((tau - t)/(sigma_t*sqrt(2)))))
  }

  A*integrate(intEq, lower=0, upper=6*t)$value
}

