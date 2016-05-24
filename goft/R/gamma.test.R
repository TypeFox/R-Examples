gammadist.test = function(x)
{  
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  n <- length(x)
  samplerange <- max(x) - min(x)
  if (samplerange == 0) 
    stop("all observations are identical")
  if( min(x) < 0 )
    stop("The dataset contains negative observations. \nAll data must be positive real numbers.")
  z <-  log(x)
  x.bar   <- mean(x)
  s2.x    <- var(x)
  b.check <- cov(x,z)  
  a.check <- x.bar/b.check
  v       <- sqrt(n*a.check)*(s2.x/(x.bar*b.check)-1)
  p.value <- 2*pnorm(abs(v),mean=0,sd=sqrt(2),lower.tail=FALSE)
  results <- list(statistic = c("V" = v), p.value = p.value, 
                  method = "Test of fit for the Gamma distribution", 
                  data.name = DNAME)
  class(results) = "htest"
  return(results)
}
