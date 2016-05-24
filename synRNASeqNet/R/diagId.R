diagId <-
function(n){
  if(n < 1) stop("n must be positive")
  if(n == 1) return(1)
  ans <- c(1, rep(0, (n-1)))
  for(i in 2:n) ans[i] <- ans[i-1] + i
  return(ans)
}
