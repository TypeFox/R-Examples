.dfmvnorm.xy <- function(x, y, w, Theta, xi, yi, ...)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    fi <- dmvnorm(cbind(x, y), mean = Theta[[i]]$theta1[c(xi, yi)], sigma = Theta[[i]]$theta2[c(xi, yi), c(xi, yi)], ...) 

    f <- f + w[i] * fi
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .dfmvnorm.xy
