propCI_agresticoull <-
function(x, n, l){

  if (x == 0){
    t <- n + qnorm(2 * l[3]) ^ 2
    p <- (x + .5 * (qnorm(2 * l[3])) ^ 2) / t
    lw <- 0
    up <- p - qnorm(2 * l[3]) * sqrt(p * (1 - p) / t)
  } else if (x == n){
    t <- n + qnorm(2 * l[3]) ^ 2
    p <- (x + .5 * (qnorm(2 * l[3])) ^ 2) / t
    lw <- p + qnorm(2 * l[3]) * sqrt(p * (1 - p) / t)
    up <- 1
  } else{
    t <- n + qnorm(l[3]) ^ 2
    p <- (x + .5 * (qnorm(l[3])) ^ 2) / t
    lw <- p + qnorm(l[3]) * sqrt(p * (1 - p) / t)
    up <- p - qnorm(l[3]) * sqrt(p * (1 - p) / t)
  }

  return(c(lw, up))
}