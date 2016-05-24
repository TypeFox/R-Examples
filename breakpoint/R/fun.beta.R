fun.beta <-
function(a, L0, L){
  x <- (a[1] - L0) / (L - L0)
  if(a[2] == 0)
    a[2] <- 0.00000001
  y <- a[2] / {(L - L0) ^ 2}
  beta.est <- (1 - x) * (((x * (1 - x)) / y) - 1)
}
