
weighted.quantile.top <- function(dat, weights, p){
  weights <- weights[order(dat)]
  n       <- length(weights)
  dat     <- sort(dat)
  q       <- sum(weights)*p
  m<- w   <-0
  for (i in n : 1){
    w <-  w+weights[i]
    if (w >= q) {m<- dat[i] ; break}
  }
  m
 }