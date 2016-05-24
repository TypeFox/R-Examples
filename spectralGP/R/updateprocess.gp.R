"updateprocess.gp" <-
function(object,...){
  # given coefficient values, calculates the process values via the inverse FFT
  # this is an internal function is not meant to be called by the user
  object$process=c(Re(fft(object$coeff,inverse=TRUE)))/sqrt(prod(object$gridsize))
  # divisor ensures that the process variance is one (it compensates
  # for scaling the coefficient variances by prod(object$gridsize) in
  # calc.variances.gp()
  return(NULL)
}
