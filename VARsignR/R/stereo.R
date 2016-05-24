stereo <-
function(v){
  llength <- length(v)
  ll <- llength + 1
  s <- matrix(NA, nrow=ll, ncol=1)
  for(i in 1:llength){
    s[i] <- 2*v[i] / (sqrt(sum(v^2))^2 +1)
  }
  s[ll] <- (sqrt(sum(v^2))^2 -1) / (sqrt(sum(v^2))^2 +1)
  return(s)
}
