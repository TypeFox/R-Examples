statVec <- function(x, mean, sd) 
{
  X <- x
  MEAN <- mean
  SD <- sd
  Z <- (((X - mean(X, na.rm = TRUE))/sd(X, na.rm = TRUE))) * SD + MEAN
  return(Z)
}


