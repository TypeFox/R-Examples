bbootstrap <-
function(xx, B = 100) {
  tmp <- rep(0, B)
  b <- 1
  while(b <= B) {
    tmp[b] <- median(sample(xx, replace = TRUE))
    b <- b + 1
  }
  return(sqrt(mean((tmp - mean(tmp))^2)))
}
