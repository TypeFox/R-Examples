.pfmvnorm.xy <- function(x, y, w, Theta, xi, yi, ...)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    for (j in 1:length(x)) {
      fi <- pmvnorm(upper = c(x[j], y[j]), mean = Theta[[i]]$theta1[c(xi, yi)], sigma = Theta[[i]]$theta2[c(xi, yi), c(xi, yi)], ...)

      f[j] <- f[j] + w[i] * fi
    }
  }  

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .pfmvnorm.xy
