thetaGT <-
function(cellCounts){
  ans <- thetaML(cellCounts)
  
  n <- sum(cellCounts)
  m1 <- sum(cellCounts == 1)
  if(m1 == n) m1 <- m1 - 1 #avoid (1 - m1/n) = 0
  
  ans <- (1 - m1/n)*ans
  return(ans)
}
