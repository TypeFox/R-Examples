df.inv <-
  function (d, df, lambda = 1, iterations = 15) 
{
  df.diriv <- function(d, lambda) -sum(d*lambda /(1 + lambda*d)^2)
  current.df <- sum(1/(1 + lambda*d))
  if (abs((df - current.df)/df) < 1e-04 | iterations == 1) 
    return(list(lambda = lambda, df = current.df))
  else {
    lambda <- exp(logb(lambda) - (current.df - df)/df.diriv(d, 
                                                            lambda))
    Recall(d, df, lambda, iterations - 1)
  }
}
