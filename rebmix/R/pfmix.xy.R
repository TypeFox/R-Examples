.pfmix.xy <- function(x, y, w, xTheta, yTheta, ...)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    if (xTheta[[i]]$pdf == .rebmix$pdf[1]) {
      fix <- pnorm(as.numeric(x), mean = as.numeric(xTheta[[i]]$theta1), sd = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[2]) {
      fix <- plnorm(as.numeric(x), meanlog = as.numeric(xTheta[[i]]$theta1), sdlog = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[3]) {
      fix <- pweibull(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[4]) {
      fix <- pbinom(as.integer(x), size = as.integer(xTheta[[i]]$theta1), prob = as.numeric(xTheta[[i]]$theta2), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[5]) {
      fix <- ppois(as.integer(x), lambda = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[6]) {
      fix <- pdirac(as.numeric(x), location = as.numeric(xTheta[[i]]$theta1), ...)
    }
    else
    if (xTheta[[i]]$pdf == .rebmix$pdf[7]) {
      fix <- pgamma(as.numeric(x), scale = as.numeric(xTheta[[i]]$theta1), shape = as.numeric(xTheta[[i]]$theta2), ...)
    }    

    if (yTheta[[i]]$pdf == .rebmix$pdf[1]) {
      fiy <- pnorm(as.numeric(y), mean = as.numeric(yTheta[[i]]$theta1), sd = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[2]) {
      fiy <- plnorm(as.numeric(y), meanlog = as.numeric(yTheta[[i]]$theta1), sdlog = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[3]) {
      fiy <- pweibull(as.numeric(y), scale = as.numeric(yTheta[[i]]$theta1), shape = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[4]) {
      fiy <- pbinom(as.integer(y), size = as.integer(yTheta[[i]]$theta1), prob = as.numeric(yTheta[[i]]$theta2), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[5]) {
      fiy <- ppois(as.integer(y), lambda = as.numeric(yTheta[[i]]$theta1), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[6]) {
      fiy <- pdirac(as.numeric(y), location = as.numeric(yTheta[[i]]$theta1), ...)
    }
    else
    if (yTheta[[i]]$pdf == .rebmix$pdf[7]) {
      fiy <- pgamma(as.numeric(y), scale = as.numeric(yTheta[[i]]$theta1), shape = as.numeric(yTheta[[i]]$theta2), ...)
    }    
    

    f <- f + w[i] * fix * fiy
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .pfmix.xy
