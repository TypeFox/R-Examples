"pchi" <-
function(q, df = 2, lower.tail = TRUE, ...)
  {
    ploum <- function(y)
      dchi(y,df)
    res <- sapply(q, function(i) integrate(ploum,0,i,...)$value)
    if (lower.tail)
      return(res)
    return(1-res)
  }

