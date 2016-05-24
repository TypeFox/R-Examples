"dchi" <-
function(x, df = 2)
  {
    r <- (2^(1-df/2))*(x^(df-1))*(exp(-(x^2)/2)) / gamma(df/2)
    return(r)
  }

