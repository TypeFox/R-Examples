LotkasC <- function(N)
{
  P <-20
  increm <- c(1:(P-1))
  sum <- sum(1/increm^N)
  part1 <- sum
  part2 <- 1/((N-1)*(P^(N-1)))
  part3 <- 1/(2*(P^N))
  part4 <- N/(24*(P-1)^(N+1))
  result <-(part1+part2+part3+part4)
  result <- 1/result
  return(result)
}
