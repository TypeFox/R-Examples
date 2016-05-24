power.set <-
function(x){
  if(length(x) == 0) return(vector(mode(x), 0))
  x <- sort(unique(x)); n <- length(x); K <- NULL
  for(m in x) K <- rbind(cbind(K, FALSE), cbind(K, TRUE))
  out <- apply(K, 1, function(x, s) s[x], s = x)
  out <- out[-c(1, length(out))]
  l <- length(out); i <- 1
  out[!sapply(out, length)>=ceiling(n/2+.5)]
}
