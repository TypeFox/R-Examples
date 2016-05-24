crossmatchdist <- function(bigN,n)
  {
  if (bigN%%2 == 1)
    {
    stop("The number of subjects, bigN, should be even")
    return(NA)
    }
  I <- bigN/2
  dist <- matrix(0,4,0)
  for (a1 in 0:I)
    {
    a2 <- (n-a1)/2
    if ( (floor(a2) == a2 ) & ( a2 >= 0 ))
      {
      a0 <- I-(a1+a2)
      if ( a0>=0 )
        {
        pr <- factorial(I)/choose(bigN,n)
        pr <- pr*( 2^a1 ) / ( factorial(a0) * factorial(a1) * factorial(a2) )
        dist <- cbind(dist, c(a0,a1,a2,pr))
        }
      }
    }
  dist <- rbind(dist,cumsum(dist[4,]))
  dist
  }