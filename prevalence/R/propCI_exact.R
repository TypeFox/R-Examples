propCI_exact <-
function(x, n, l){

  if (x == 0){
    lw <- 0
    up <- 1 - l[3] ^ (1/n)
  } else if (x == n){
    lw <- l[3] ^ (1/n)
    up <- 1
  } else{
    lw <- qbeta(l[1], x, n - x + 1, lower.tail = T)
    up <- qbeta(l[2], x + 1, n - x, lower.tail = T)
  }

  return(c(lw, up))
}