initialise  <- function(x,n)
{
  sam <- sample(seq(1:nrow(x)),n)
  xpop <- x[sam,]
  return(xpop)
}

