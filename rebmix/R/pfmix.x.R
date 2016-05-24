.pfmix.x <- function(x, w, xTheta, ...)
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

    f <- f + w[i] * fix
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .pfmix.x
