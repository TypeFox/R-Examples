rmnorm <-
function (n = 1,
                    mean = rep(0, d),
                    varcov)
{
  d <- if (is.matrix(varcov))
    ncol(varcov)
  else
    1
  #Handle non-positive definite 'varcov': Due to numerical rounding issues
  #there is sometimes a problem with varcov not beeing pos. def. This is
  #due to eigenvalues near 0 and with some numerical rounding very small
  #(~10^20) negative eigenvalues might appear.
  z <- try(matrix(rnorm(n * d), n, d) %*% chol(varcov), silent = TRUE)
  if(inherits(z, "try-error"))
    z <- try(matrix(rnorm(n * d), n, d) %*% chol(nearPD(varcov)), silent = TRUE)
  if(inherits(z, "try-error"))
    return(rep(NA, d))
  y <- t(mean + t(z))
  if(n == 1 & !is.null(colnames(varcov)))
  {
    y <- c(y)
    names(y) <- colnames(varcov)
  }
  return(y)
}
