ig.test <- function(x)
{  
  DNAME <- deparse(substitute(x))
  stopifnot(is.numeric(x))
  samplerange <- max(x) - min(x)
  if (samplerange == 0) 
    stop("all observations are identical")
  if( min(x) < 0 )
    stop("The dataset contains negative observations. \nAll data must be non-negative real numbers.")
  z <- ((x - mean(x))**2)/x  # transformation given in equation (6)
  res <- .gammadist.test2(z)
  ad <- res$statistic
  p.value <- res$p.value
    results <- list(statistic =  ad, p.value = p.value, method = "A transformation test for Inverse Gaussian distributions", 
                  data.name = DNAME)
  class(results) = "htest"
  return(results)
}