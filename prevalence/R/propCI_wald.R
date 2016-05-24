propCI_wald <-
function(x, n, l){
  p <- x / n

  if (x == 0){
    lw <- 0
    up <- p - qnorm(l[3]) * sqrt(p * (1 - p) / n)
  } else if (x == n){
    lw <- p + qnorm(l[3]) * sqrt(p * (1 - p) / n)
    up <- 1
  } else{
    lw <- p + qnorm(l[3]) * sqrt(p * (1 - p) / n)
    up <- p - qnorm(l[3]) * sqrt(p * (1 - p) / n)
  }

  return(c(lw, up))
}
