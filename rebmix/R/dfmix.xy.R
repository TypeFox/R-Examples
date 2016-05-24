.dfmix.xy <- function(x, y, w, xTheta, yTheta, ...)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    if (xTheta[[i]]$pdf == .rebmix$pdf[1]) {
      fix <- dnorm(as.numeric(x), mean = as.numeric(xTheta[[i]]$theta1), sd = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[2]) {
      fix <- dlnorm(as.numeric(x), meanlog = as.numeric(xTheta[[i]]$theta1), sdlog = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[3]) {
      fix <- dweibull(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[4]) {
      fix <- dbinom(as.integer(x), size = as.integer(xTheta[[i]]$theta1), prob = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[5]) {
      fix <- dpois(as.integer(x), lambda = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[6]) {
      fix <- ddirac(as.numeric(x), location = as.numeric(xTheta[[i]]$theta1))
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[7]) {
      fix <- dgamma(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }    

    if (yTheta[[i]]$pdf == .rebmix$pdf[1]) {
      fiy <- dnorm(as.numeric(y), mean = as.numeric(yTheta[[i]]$theta1), sd = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[2]) {
      fiy <- dlnorm(as.numeric(y), meanlog = as.numeric(yTheta[[i]]$theta1), sdlog = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[3]) {
      fiy <- dweibull(as.numeric(y), scale = as.numeric(yTheta[[i]]$theta1), shape = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[4]) {
      fiy <- dbinom(as.integer(y), size = as.integer(yTheta[[i]]$theta1), prob = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[5]) {
      fiy <- dpois(as.integer(y), lambda = as.numeric(yTheta[[i]]$theta1), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[6]) {
      fiy <- ddirac(as.numeric(y), location = as.numeric(yTheta[[i]]$theta1))
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[7]) {
      fiy <- dgamma(as.numeric(y), scale = as.numeric(yTheta[[i]]$theta1), shape = as.numeric(yTheta[[i]]$theta2), ...)
    }    

    f <- f + w[i] * fix * fiy
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .dfmix.xy
