lineslope <-
function(K) {
  line.slope=c()
  for (i in 1:K)
  {
    line.slope[i]=-tan((K+1-i)*pi/(2*(K+1)))
  }
  return(line.slope)
}
