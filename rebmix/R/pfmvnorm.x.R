.pfmvnorm.x <- function(x, w, Theta, xi, ...)
{
  f <- array(data = 0.0, dim = length(x), dimnames = NULL)

  for (i in 1:length(w)) {
    fi <- pnorm(as.numeric(x), mean = as.numeric(Theta[[i]]$theta1[xi]), sd = as.numeric(Theta[[i]]$theta2[xi]^0.5), ...)

    f <- f + w[i] * fi
  }

  rm(list = ls()[!(ls() %in% c("f"))])

  return(f)
} ## .pfmvnorm.x
