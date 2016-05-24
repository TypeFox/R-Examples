dgenpois <-
function(x, lambda1, lambda2)
{ 
	
  if (length(x) < max(length(lambda1), length(lambda2)))
  {
  	x <- c(rep(x, times = max(length(lambda1), length(lambda2))))
  }
  if (length(lambda1) < max(length(x), length(lambda2)))
  {
  	lambda1 <- c(rep(lambda1, times = max(length(x), length(lambda2))))
  }
  if (length(lambda2) < max(length(x), length(lambda1)))
  {
  	lambda2 <- c(rep(lambda2, times = max(length(x), length(lambda1))))
  }
  a <- NULL
  for (j in 1:max(length(x), length(lambda1), length(lambda2)))
  { 
  	if(x[j] < 2)
  { 
  	b = (lambda1[j] * (lambda1[j] + x[j] * lambda2[j])^(x[j] - 1) * (exp(-(
    lambda1[j] + x[j] * lambda2[j])))) / (factorial(x[j]))
  } else  { 
  	  f1 <- (lambda1[j] + x[j] * lambda2[j])
      g1 <- exp(-lambda2[j])
      e <- 1
      for (i in 2:x[j])
      { 
      	e <- e * ((f1 * g1) / i)
      }
      d1 <- lambda1[j] * exp(-lambda1[j]) * g1
      b <- e * d1
      if (is.na(b))
      {
      	b <- 4.940656e-324
      }
      if(b == Inf)
      {
      	b <- 4.940656e-324
      }
      if(b == -Inf)
      {
      	b <- 4.940656e-324
      }
    }
    a <- c(a,b)
  }
return(a)
}
