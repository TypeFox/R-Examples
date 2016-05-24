"Normalize" <- 
function(x)
{

uniformize <- function (x) 
{ 
	x <- rank(x, 
             na.last = "keep", 
             ties.method = "average")
    n <- sum(!is.na(x))
    x / (n + 1)
  }
  
  return( qnorm(uniformize(x), mean=mean(x, na.rm=TRUE), sd=sd(x, na.rm=TRUE) ) )
  }
  
  
 