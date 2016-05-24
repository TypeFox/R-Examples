decToBin <-
function(x, m) {

  bin <- array(0, c(length(x),m))
  for(i in 1:m) {
    bin[,m - i + 1] <- x %% 2
    x <- x %/% 2
  }

  return(bin)

}
