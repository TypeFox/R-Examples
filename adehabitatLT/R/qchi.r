"qchi" <-
function(p,df = 2, lower.tail = TRUE)
  {
    return(sqrt(qchisq(p, df = df, lower.tail=lower.tail)))
  }

