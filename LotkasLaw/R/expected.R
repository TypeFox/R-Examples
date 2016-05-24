expected <- function(Table,C,N)
{
  value <- Table[,1:1]^N
  nvalue <- 1/value
  part2 <- C*nvalue
  return(part2)
  
}
