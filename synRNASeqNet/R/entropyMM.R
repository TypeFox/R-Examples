entropyMM <-
function(cellCounts, unit = unit){
  ans <- entropyML(cellCounts, unit = unit)
  
  n <- sum(cellCounts)
  m <- sum(cellCounts > 0)
  
  ans <- ans + (m - 1)/(2*n)
  return(ans)
}
