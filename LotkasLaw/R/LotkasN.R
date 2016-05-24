LotkasN <- function(Sums,FullTable)
{
  N <- nrow(FullTable)
  lx <- Sums[3]
  ly <- Sums[4]
  xy <- Sums[5]
  x2 <- Sums[6]
  lx2 <- lx^2
  top <- (N*xy) - (lx*ly)
  bottom <- (N*x2) - (lx2)
  Nfinal <- top/bottom
  return(Nfinal)
}
