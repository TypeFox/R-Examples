qfun <-
  function(x){# define Q function
    x=exp(.5*x^2+pnorm(x,lower.tail = FALSE, log.p = TRUE))
  }
